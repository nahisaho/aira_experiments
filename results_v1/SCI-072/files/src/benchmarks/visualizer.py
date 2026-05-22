from __future__ import annotations

import sys
sys.path.insert(0, '/home/nahisaho/GitHub/aira/projects/f78a410a-e891-4bfe-88c5-74583c82bef3/workspace')

from pathlib import Path
from typing import Any, Sequence

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from src.core.solution import Path as SolverPath
from src.core.solution import Solution as SolverSolution

from .types import GridEnvironment, Position


class Visualizer:
    """Matplotlib visualizations for benchmark outcomes."""

    def __init__(self, figures_dir: str = 'figures') -> None:
        self.figures_dir = Path(figures_dir)
        self.figures_dir.mkdir(parents=True, exist_ok=True)

    def plot_scalability(self, results_df: pd.DataFrame, metric: str = 'runtime', save_path: str | None = None) -> Path:
        """Plot runtime or another metric against the number of agents."""
        metric_column = 'runtime_seconds' if metric == 'runtime' else metric
        if metric_column not in results_df.columns or 'num_agents' not in results_df.columns:
            raise ValueError('Scalability plotting requires num_agents and the requested metric column.')

        fig, ax = plt.subplots(figsize=(8, 5))
        solvers = sorted(results_df['solver'].unique()) if 'solver' in results_df.columns else ['solver']
        colors = plt.cm.viridis(np.linspace(0.15, 0.85, len(solvers)))
        for color, solver_name in zip(colors, solvers):
            subset = results_df if solver_name == 'solver' else results_df[results_df['solver'] == solver_name]
            subset = subset.sort_values('num_agents')
            ax.plot(subset['num_agents'], subset[metric_column], marker='o', color=color, label=solver_name)
        ax.set_xlabel('Number of agents')
        ax.set_ylabel(metric_column.replace('_', ' ').title())
        ax.set_title('Scalability Benchmark')
        if 'solver' in results_df.columns:
            ax.legend()
        ax.grid(True, alpha=0.3)
        return self._finalize_figure(fig, save_path, 'scalability.png')

    def plot_quality_comparison(self, results_df: pd.DataFrame, save_path: str | None = None) -> Path:
        """Plot bar chart of mean sum-of-cost ratios."""
        ratio_column = 'suboptimality_ratio' if 'suboptimality_ratio' in results_df.columns else 'soc_ratio'
        if ratio_column not in results_df.columns or 'solver' not in results_df.columns:
            raise ValueError('Quality comparison requires solver and suboptimality ratio columns.')

        grouped = results_df.groupby('solver', dropna=False)[ratio_column].mean().reset_index()
        fig, ax = plt.subplots(figsize=(8, 5))
        colors = plt.cm.cividis(np.linspace(0.2, 0.85, len(grouped)))
        ax.bar(grouped['solver'], grouped[ratio_column], color=colors)
        ax.set_xlabel('Solver')
        ax.set_ylabel('Mean SOC Ratio')
        ax.set_title('Quality Comparison')
        ax.grid(axis='y', alpha=0.3)
        return self._finalize_figure(fig, save_path, 'quality_comparison.png')

    def plot_solution(self, env: GridEnvironment, solution: Any, save_path: str | None = None) -> Path:
        """Plot a static path overlay on the environment grid."""
        paths = self._extract_paths(solution)
        fig, ax = plt.subplots(figsize=(8, 8))
        grid = self._environment_grid(env)
        ax.imshow(grid, cmap='Greys', origin='upper')
        colors = plt.cm.viridis(np.linspace(0.1, 0.95, max(len(paths), 1)))
        for index, path in enumerate(paths):
            if not path:
                continue
            xs = [cell[0] for cell in path]
            ys = [cell[1] for cell in path]
            ax.plot(xs, ys, color=colors[index], linewidth=2, alpha=0.9)
            ax.scatter(xs[0], ys[0], color=colors[index], marker='o', s=40)
            ax.scatter(xs[-1], ys[-1], color=colors[index], marker='x', s=50)
        if env.stations:
            station_x, station_y = zip(*env.stations)
            ax.scatter(station_x, station_y, color='gold', marker='s', s=45, label='Stations')
            ax.legend()
        ax.set_xlabel('X coordinate')
        ax.set_ylabel('Y coordinate')
        ax.set_title('Solution Paths')
        ax.set_aspect('equal')
        return self._finalize_figure(fig, save_path, 'solution.png')

    def plot_warehouse_heatmap(self, env: GridEnvironment, solution: Any, save_path: str | None = None) -> Path:
        """Plot cell visit counts for a warehouse solution."""
        heatmap = np.zeros((env.height, env.width), dtype=float)
        for path in self._extract_paths(solution):
            for x, y in path:
                heatmap[y, x] += 1.0
        for x, y in env.obstacles:
            heatmap[y, x] = np.nan
        fig, ax = plt.subplots(figsize=(8, 6))
        image = ax.imshow(heatmap, cmap='cividis', origin='upper')
        ax.set_xlabel('X coordinate')
        ax.set_ylabel('Y coordinate')
        ax.set_title('Warehouse Congestion Heatmap')
        fig.colorbar(image, ax=ax, label='Visits')
        return self._finalize_figure(fig, save_path, 'warehouse_heatmap.png')

    def plot_convergence(self, results_df: pd.DataFrame, save_path: str | None = None) -> Path:
        """Plot convergence over iteration or step."""
        x_column = 'iteration' if 'iteration' in results_df.columns else 'step'
        y_column = 'objective' if 'objective' in results_df.columns else 'cost'
        if x_column not in results_df.columns or y_column not in results_df.columns:
            raise ValueError('Convergence plotting requires iteration/step and objective/cost columns.')

        fig, ax = plt.subplots(figsize=(8, 5))
        solvers = sorted(results_df['solver'].unique()) if 'solver' in results_df.columns else ['solver']
        colors = plt.cm.viridis(np.linspace(0.15, 0.85, len(solvers)))
        for color, solver_name in zip(colors, solvers):
            subset = results_df if solver_name == 'solver' else results_df[results_df['solver'] == solver_name]
            subset = subset.sort_values(x_column)
            ax.plot(subset[x_column], subset[y_column], color=color, label=solver_name)
        ax.set_xlabel(x_column.replace('_', ' ').title())
        ax.set_ylabel(y_column.replace('_', ' ').title())
        ax.set_title('Solver Convergence')
        if 'solver' in results_df.columns:
            ax.legend()
        ax.grid(True, alpha=0.3)
        return self._finalize_figure(fig, save_path, 'convergence.png')

    def plot_lifelong_throughput(self, results_df: pd.DataFrame, save_path: str | None = None) -> Path:
        """Plot lifelong throughput over time or by solver."""
        fig, ax = plt.subplots(figsize=(8, 5))
        if {'time', 'throughput'} <= set(results_df.columns):
            solvers = sorted(results_df['solver'].unique()) if 'solver' in results_df.columns else ['solver']
            colors = plt.cm.cividis(np.linspace(0.15, 0.85, len(solvers)))
            for color, solver_name in zip(colors, solvers):
                subset = results_df if solver_name == 'solver' else results_df[results_df['solver'] == solver_name]
                subset = subset.sort_values('time')
                ax.plot(subset['time'], subset['throughput'], marker='o', color=color, label=solver_name)
            if 'solver' in results_df.columns:
                ax.legend()
            ax.set_xlabel('Time')
        else:
            if 'throughput' not in results_df.columns or 'solver' not in results_df.columns:
                raise ValueError('Throughput plotting requires either time+throughput or solver+throughput columns.')
            grouped = results_df.groupby('solver', dropna=False)['throughput'].mean().reset_index()
            colors = plt.cm.cividis(np.linspace(0.2, 0.85, len(grouped)))
            ax.bar(grouped['solver'], grouped['throughput'], color=colors)
            ax.set_xlabel('Solver')
        ax.set_ylabel('Throughput')
        ax.set_title('Lifelong Throughput')
        ax.grid(True, alpha=0.3)
        return self._finalize_figure(fig, save_path, 'lifelong_throughput.png')

    def _extract_paths(self, solution: Any) -> list[list[Position]]:
        if solution is None:
            return []
        if isinstance(solution, SolverSolution):
            return [self._normalize_path(path.states) for path in solution.paths.values()]
        if isinstance(solution, dict):
            if 'solution' in solution:
                return self._extract_paths(solution['solution'])
            if 'paths' in solution:
                solution = solution['paths']
            if isinstance(solution, dict):
                return [self._normalize_path(path) for path in solution.values()]
        return [self._normalize_path(path) for path in solution]

    def _normalize_path(self, path: Sequence[Any] | SolverPath) -> list[Position]:
        if isinstance(path, SolverPath):
            return [tuple(step) for step in path.states]
        normalized: list[Position] = []
        for cell in path:
            if isinstance(cell, dict):
                normalized.append((int(cell['x']), int(cell['y'])))
            else:
                x, y = cell
                normalized.append((int(x), int(y)))
        return normalized

    def _environment_grid(self, env: GridEnvironment) -> np.ndarray:
        grid = np.zeros((env.height, env.width), dtype=float)
        for x, y in env.obstacles:
            grid[y, x] = 1.0
        return grid

    def _finalize_figure(self, fig: plt.Figure, save_path: str | None, default_name: str) -> Path:
        path = Path(save_path) if save_path else self.figures_dir / default_name
        path.parent.mkdir(parents=True, exist_ok=True)
        fig.tight_layout()
        fig.savefig(path, dpi=300, bbox_inches='tight')
        plt.close(fig)
        return path
