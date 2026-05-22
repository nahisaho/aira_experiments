"""
ROS Launch Configuration (Python-based for ROS2)
Launches all nodes for the semi-autonomous suturing system.
"""

# ROS2 launch file equivalent (for documentation purposes)
# In production, this would be a proper ROS2 launch.py

LAUNCH_CONFIG = {
    'nodes': [
        {
            'package': 'suturing_system',
            'executable': 'lfd_node',
            'name': 'lfd_trajectory_generator',
            'parameters': {
                'method': 'gmr',
                'n_components': 5,
                'demo_directory': 'data/demonstrations/',
            },
            'remappings': [
                ('trajectory_out', '/suturing/reference_trajectory'),
            ],
        },
        {
            'package': 'suturing_system',
            'executable': 'tissue_model_node',
            'name': 'tissue_deformation_model',
            'parameters': {
                'model_type': 'mass_spring',
                'update_rate': 1000,  # Hz
            },
            'remappings': [
                ('deformation_out', '/suturing/tissue_deformation'),
                ('force_in', '/suturing/tool_force'),
            ],
        },
        {
            'package': 'suturing_system',
            'executable': 'compliance_node',
            'name': 'compliance_controller',
            'parameters': {
                'mode': 'impedance',
                'update_rate': 1000,
            },
            'remappings': [
                ('force_in', '/suturing/measured_force'),
                ('pose_cmd_out', '/suturing/compliant_pose'),
            ],
        },
        {
            'package': 'suturing_system',
            'executable': 'visual_servo_node',
            'name': 'visual_servo_controller',
            'parameters': {
                'mode': 'pbvs',
                'gain': 0.5,
            },
            'remappings': [
                ('image_left', '/endoscope/left/image_raw'),
                ('image_right', '/endoscope/right/image_raw'),
                ('velocity_cmd', '/suturing/vs_velocity'),
            ],
        },
        {
            'package': 'suturing_system',
            'executable': 'safety_node',
            'name': 'safety_monitor',
            'parameters': {
                'force_limit_critical': 10.0,
                'workspace_radius': 0.15,
                'update_rate': 1000,
            },
            'remappings': [
                ('safety_status', '/suturing/safety_state'),
                ('emergency_stop', '/dvrk/emergency_stop'),
            ],
        },
        {
            'package': 'suturing_system',
            'executable': 'coordinator_node',
            'name': 'suturing_coordinator',
            'parameters': {
                'config_file': 'config/suturing_config.yaml',
            },
            'remappings': [
                ('phase_status', '/suturing/phase_status'),
            ],
        },
    ],
    'topics': {
        '/suturing/reference_trajectory': 'geometry_msgs/PoseArray',
        '/suturing/tissue_deformation': 'sensor_msgs/PointCloud2',
        '/suturing/measured_force': 'geometry_msgs/WrenchStamped',
        '/suturing/compliant_pose': 'geometry_msgs/PoseStamped',
        '/suturing/vs_velocity': 'geometry_msgs/TwistStamped',
        '/suturing/safety_state': 'std_msgs/String',
        '/suturing/phase_status': 'std_msgs/String',
        '/endoscope/left/image_raw': 'sensor_msgs/Image',
        '/endoscope/right/image_raw': 'sensor_msgs/Image',
        '/dvrk/PSM1/position_cartesian_current': 'geometry_msgs/PoseStamped',
        '/dvrk/PSM1/state_joint_current': 'sensor_msgs/JointState',
    },
}


def print_launch_info():
    """Print launch configuration summary."""
    print("=" * 60)
    print("Suturing System ROS2 Launch Configuration")
    print("=" * 60)
    for node in LAUNCH_CONFIG['nodes']:
        print(f"\n  Node: {node['name']}")
        print(f"    Package: {node['package']}")
        print(f"    Executable: {node['executable']}")
        if 'parameters' in node:
            for k, v in node['parameters'].items():
                print(f"    Param {k}: {v}")

    print(f"\n  Topics ({len(LAUNCH_CONFIG['topics'])}):")
    for topic, msg_type in LAUNCH_CONFIG['topics'].items():
        print(f"    {topic} [{msg_type}]")


if __name__ == '__main__':
    print_launch_info()
