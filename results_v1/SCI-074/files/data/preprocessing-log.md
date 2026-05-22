# Preprocessing Log

- Input data source assumed: Intel RealSense D455 stereo infrared images, D455 depth point clouds, BMI088 IMU.
- No raw mission data was processed in this task; this file documents design-time preprocessing assumptions.
- Assumed preprocessing chain for runtime deployment:
  1. camera rectification from offline calibration,
  2. timestamp normalization into a common monotonic clock,
  3. depth validity filtering and point-cloud cropping,
  4. optional statistical outlier removal before map fusion,
  5. Allan variance-based IMU noise identification for final `vio_config.yaml` tuning.
