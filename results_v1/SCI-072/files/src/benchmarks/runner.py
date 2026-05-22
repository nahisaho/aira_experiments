from __future__ import annotations

import sys
sys.path.insert(0, '/home/nahisaho/GitHub/aira/projects/f78a410a-e891-4bfe-88c5-74583c82bef3/workspace')

import copy
from contextlib import contextmanager
import inspect
from pathlib import Path
import signal
import threading
import tracemalloc
from typing import Any, Callable, Iterable, Sequence

import pandas as pd

from src.utils import Timer, setup_logger
from src.utils.logger import ProcessLogger

from .metrics import MetricsCollector
from .scenarios import ScenarioGenerator
from .types import Agent, GridEnvironment


class BenchmarkRunner:
    """Orchestrate solver benchmark experiments."""

    def __init__(self, output_dir: str = 'results') -> None:
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.solvers: list[dict[str, Any]] = []
        self.scenarios: list[dict[str, Any]] = []
        self.logger = setup_logger(self.__class__.__name__)
        self.process_logger = ProcessLogger()

    def add_solver(self, name: str, solver_class: type | Callable[..., Any], **kwargs: Any) -> None:
        """Register a solver implementation."""
        self.solvers.append({'name': name, 'solver_class': solver_class, 'kwargs': dict(kwargs)})

    def add_scenario(self, name: str, env: GridEnvironment, agents: Sequence[Agent], **kwargs: Any) -> None:
        """Register a benchmark scenario."""
        self.scenarios.append({'name': name, 'env': env, 'agents': list(agents), 'kwargs': dict(kwargs)})

    def run_all(self, repeats: int = 3) -> pd.DataFrame:
        """Run every registered solver on every registered scenario."""
        records: list[dict[str, Any]] = []
        for solver_spec in self.solvers:
            for scenario_spec in self.scenarios:
                for repeat in range(1, repeats + 1):
                    records.append(self._execute_run(solver_spec, scenario_spec, repeat=repeat))
        return pd.DataFrame(records)

    def run_scalability_test(
        self,
        solver_class: type | Callable[..., Any],
        env_generator: Callable[[int], Any],
        agent_counts: list[int] | None = None,
        timeout: int = 300,
    ) -> pd.DataFrame:
        """Run one solver across a range of agent counts."""
        counts = agent_counts or [10, 20, 50, 100, 200, 500, 1000]
        solver_spec = self._normalize_solver_entry(solver_class)
        records: list[dict[str, Any]] = []
        for count in counts:
            generated = env_generator(count)
            env, agents, scenario_kwargs = self._normalize_generated_case(generated, count)
            scenario_spec = {
                'name': f'scalability_{count}',
                'env': env,
                'agents': agents,
                'kwargs': scenario_kwargs,
            }
            record = self._execute_run(solver_spec, scenario_spec, repeat=1, timeout=timeout)
            record['num_agents'] = count
            records.append(record)
        return pd.DataFrame(records)

    def run_quality_comparison(
        self,
        optimal_solver: Any,
        approx_solvers: Any,
        scenarios: Sequence[Any],
    ) -> pd.DataFrame:
        """Compare approximate solvers against an optimal solver."""
        optimal_spec = self._normalize_solver_entry(optimal_solver, default_name='optimal')
        approx_specs = [self._normalize_solver_entry(entry) for entry in self._iter_solver_entries(approx_solvers)]
        scenario_specs = [self._normalize_scenario_entry(entry) for entry in scenarios]

        records: list[dict[str, Any]] = []
        for scenario_spec in scenario_specs:
            optimal_record = self._execute_run(optimal_spec, scenario_spec, repeat=1)
            optimal_cost = optimal_record.get('sum_of_costs')
            optimal_record['comparison_type'] = 'optimal'
            records.append(optimal_record)
            for solver_spec in approx_specs:
                record = self._execute_run(solver_spec, scenario_spec, repeat=1)
                record['optimal_cost'] = optimal_cost
                record['suboptimality_ratio'] = (
                    record['sum_of_costs'] / optimal_cost if optimal_cost not in (None, 0) else None
                )
                record['comparison_type'] = 'approximate'
                records.append(record)
        return pd.DataFrame(records)

    def save_results(self, df: pd.DataFrame, filename: str) -> Path:
        """Save results to a CSV file."""
        path = self.output_dir / filename
        path.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(path, index=False)
        return path

    def generate_summary(self, df: pd.DataFrame) -> str:
        """Create a markdown table summarizing benchmark results."""
        if df.empty:
            return 'No benchmark results available.'

        group_columns = [column for column in ('solver', 'scenario') if column in df.columns]
        if not group_columns:
            raise ValueError('Results DataFrame must contain at least a solver or scenario column.')

        summary = (
            df.groupby(group_columns, dropna=False)
            .agg(
                runs=('success', 'size') if 'success' in df.columns else (df.columns[0], 'size'),
                success_rate=('success', 'mean') if 'success' in df.columns else (df.columns[0], 'size'),
                mean_runtime_seconds=('runtime_seconds', 'mean') if 'runtime_seconds' in df.columns else (df.columns[0], 'size'),
                mean_sum_of_costs=('sum_of_costs', 'mean') if 'sum_of_costs' in df.columns else (df.columns[0], 'size'),
                mean_makespan=('makespan', 'mean') if 'makespan' in df.columns else (df.columns[0], 'size'),
            )
            .reset_index()
        )
        if 'success_rate' in summary.columns:
            summary['success_rate'] = (summary['success_rate'] * 100).round(1)
        for column in summary.columns:
            if summary[column].dtype.kind in {'f', 'i'}:
                summary[column] = summary[column].map(lambda value: round(float(value), 3))

        headers = list(summary.columns)
        lines = [
            '| ' + ' | '.join(headers) + ' |',
            '| ' + ' | '.join(['---'] * len(headers)) + ' |',
        ]
        for _, row in summary.iterrows():
            lines.append('| ' + ' | '.join(str(row[column]) for column in headers) + ' |')
        return '\\n'.join(lines)

    def _execute_run(
        self,
        solver_spec: dict[str, Any],
        scenario_spec: dict[str, Any],
        repeat: int,
        timeout: int | None = None,
    ) -> dict[str, Any]:
        env = copy.deepcopy(scenario_spec['env'])
        agents = copy.deepcopy(scenario_spec['agents'])
        scenario_kwargs = copy.deepcopy(scenario_spec.get('kwargs', {}))
        solver = self._instantiate_solver(solver_spec, env, agents, timeout=timeout)

        self.process_logger.log(
            phase='execute',
            event_type='handoff_started',
            skill_or_tool='BenchmarkRunner',
            handoff_in={'solver': solver_spec['name'], 'scenario': scenario_spec['name'], 'repeat': repeat},
        )

        solution: Any = None
        raw_stats: dict[str, Any] = {}
        error_message: str | None = None
        tracemalloc.start()
        with Timer() as timer:
            try:
                with self._timeout_guard(timeout):
                    result = self._call_solver(solver, env, agents, **scenario_kwargs)
                    solution, raw_stats = self._unpack_solver_result(result, solver)
            except TimeoutError:
                raw_stats = {'success': False, 'timeout': True}
                error_message = f'Timed out after {timeout} seconds.'
            except Exception as exc:  # pragma: no cover
                raw_stats = {'success': False, 'error': str(exc)}
                error_message = str(exc)
        _, peak_memory = tracemalloc.get_traced_memory()
        tracemalloc.stop()

        raw_stats['runtime_seconds'] = raw_stats.get('runtime_seconds', timer.elapsed)
        raw_stats['memory_usage_mb'] = max(float(raw_stats.get('memory_usage_mb', 0.0) or 0.0), peak_memory / (1024 * 1024))
        metrics = MetricsCollector.compute_all(solution, env, agents, raw_stats)
        record = {
            'solver': solver_spec['name'],
            'scenario': scenario_spec['name'],
            'repeat': repeat,
            'num_agents': len(agents),
            **metrics,
        }
        if error_message:
            record['error'] = error_message
            self.logger.warning('%s on %s failed: %s', solver_spec['name'], scenario_spec['name'], error_message)
        else:
            self.logger.info('Completed %s on %s (repeat %s).', solver_spec['name'], scenario_spec['name'], repeat)

        self.process_logger.log(
            phase='execute',
            event_type='handoff_completed',
            skill_or_tool='BenchmarkRunner',
            handoff_out={'solver': solver_spec['name'], 'scenario': scenario_spec['name'], 'repeat': repeat, 'success': record['success']},
            status='ok' if record['success'] else 'error',
        )
        return record

    def _instantiate_solver(
        self,
        solver_spec: dict[str, Any],
        env: GridEnvironment,
        agents: Sequence[Agent],
        timeout: int | None = None,
    ) -> Any:
        solver_factory = solver_spec['solver_class']
        kwargs = dict(solver_spec.get('kwargs', {}))
        if timeout is not None:
            kwargs.setdefault('timeout', timeout)

        if inspect.isclass(solver_factory):
            try:
                return solver_factory(env, list(agents), **kwargs)
            except TypeError:
                return solver_factory(**kwargs)
        return solver_factory

    def _call_solver(self, solver: Any, env: GridEnvironment, agents: Sequence[Agent], **kwargs: Any) -> Any:
        if hasattr(solver, 'solve'):
            try:
                return solver.solve(env, list(agents), **kwargs)
            except TypeError:
                return solver.solve()
        if callable(solver):
            try:
                return solver(env, list(agents), **kwargs)
            except TypeError:
                return solver()
        raise TypeError('Solver must be callable or expose a solve() method.')

    def _unpack_solver_result(self, result: Any, solver: Any) -> tuple[Any, dict[str, Any]]:
        if isinstance(result, tuple):
            if len(result) == 2:
                solution, stats = result
                return solution, dict(stats or {})
            if len(result) == 1:
                return result[0], {}
        if isinstance(result, dict) and 'solution' in result:
            return result['solution'], dict(result.get('stats', {}))
        return result, dict(getattr(solver, 'stats', {}) or {})

    def _normalize_generated_case(self, generated: Any, num_agents: int) -> tuple[GridEnvironment, list[Agent], dict[str, Any]]:
        if isinstance(generated, dict):
            env = generated['env']
            agents = generated.get('agents') or ScenarioGenerator.generate_random(env, num_agents)
            kwargs = generated.get('kwargs', {})
            return env, list(agents), dict(kwargs)
        if isinstance(generated, tuple):
            if len(generated) == 3:
                env, agents, kwargs = generated
                return env, list(agents), dict(kwargs)
            if len(generated) == 2:
                env, agents = generated
                return env, list(agents), {}
        if isinstance(generated, GridEnvironment):
            return generated, ScenarioGenerator.generate_random(generated, num_agents), {}
        raise TypeError('env_generator must return GridEnvironment, (env, agents), or (env, agents, kwargs).')

    def _normalize_solver_entry(self, entry: Any, default_name: str | None = None) -> dict[str, Any]:
        if isinstance(entry, dict) and 'solver_class' in entry:
            return {
                'name': entry.get('name', default_name or 'solver'),
                'solver_class': entry['solver_class'],
                'kwargs': dict(entry.get('kwargs', {})),
            }
        if isinstance(entry, tuple):
            if len(entry) == 3:
                name, solver_class, kwargs = entry
                return {'name': name, 'solver_class': solver_class, 'kwargs': dict(kwargs)}
            if len(entry) == 2:
                name, solver_class = entry
                return {'name': name, 'solver_class': solver_class, 'kwargs': {}}
        name = default_name or getattr(entry, '__name__', entry.__class__.__name__)
        return {'name': name, 'solver_class': entry, 'kwargs': {}}

    def _normalize_scenario_entry(self, entry: Any) -> dict[str, Any]:
        if isinstance(entry, dict) and {'name', 'env', 'agents'} <= set(entry.keys()):
            return {
                'name': entry['name'],
                'env': entry['env'],
                'agents': list(entry['agents']),
                'kwargs': dict(entry.get('kwargs', {})),
            }
        if isinstance(entry, tuple):
            if len(entry) == 4:
                name, env, agents, kwargs = entry
                return {'name': name, 'env': env, 'agents': list(agents), 'kwargs': dict(kwargs)}
            if len(entry) == 3:
                name, env, agents = entry
                return {'name': name, 'env': env, 'agents': list(agents), 'kwargs': {}}
        raise TypeError('Scenarios must be dicts or tuples of (name, env, agents[, kwargs]).')

    def _iter_solver_entries(self, approx_solvers: Any) -> list[Any]:
        if isinstance(approx_solvers, dict):
            return [(name, solver) for name, solver in approx_solvers.items()]
        return list(approx_solvers)

    @contextmanager
    def _timeout_guard(self, seconds: int | None):
        if not seconds or seconds <= 0 or not hasattr(signal, 'SIGALRM') or threading.current_thread() is not threading.main_thread():
            yield
            return

        def handler(signum: int, frame: Any) -> None:
            raise TimeoutError('Solver execution exceeded time limit.')

        previous_handler = signal.signal(signal.SIGALRM, handler)
        signal.setitimer(signal.ITIMER_REAL, float(seconds))
        try:
            yield
        finally:
            signal.setitimer(signal.ITIMER_REAL, 0)
            signal.signal(signal.SIGALRM, previous_handler)
