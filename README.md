# SAR Preprocessing Toolkit

A Python toolkit for Synthetic Aperture Radar (SAR) image preprocessing, providing radiometric calibration, geometric correction, and orthorectification capabilities.

## Overview

This toolkit is designed for preprocessing SAR satellite imagery, with primary support for LT-1 satellite data. It provides a complete pipeline from raw SAR data to geocoded, radiometrically calibrated products.

## Features

- **Radiometric Calibration**: Convert raw digital numbers (DN) to calibrated backscatter coefficients (σ₀)
- **Geometric Correction**: Transform slant-range imagery to ground-range using RPC (Rational Polynomial Coefficients) models
- **Orthorectification**: Generate orthorectified products with customizable spatial resolution and projection
- **Polarimetric Processing**: Support for single and dual-polarization SAR data
- **Incidence Angle Calculation**: Compute local incidence angles for terrain correction
- **Tiled Processing**: Memory-efficient block-based processing for large images

## Supported Satellites

- **LT-1** (Luojia-1): Chinese commercial SAR satellite constellation
- Extensible architecture for other SAR sensors

## Installation

### Requirements

- Python 3.7+
- GDAL/OGR libraries

### Dependencies

```bash
pip install -r requirements.txt
```

Core dependencies:
- `numpy` - Numerical computing
- `scipy` - Scientific computing
- `rasterio` - Geospatial raster I/O
- `h5py` - HDF5 file handling
- `matplotlib` - Visualization
- `opencv-python` - Image processing
- `pyproj` - Cartographic projections
- `scikit-image` - Image processing algorithms
- `GDAL` - Geospatial data abstraction library
- `pillow` - Image file handling

## Quick Start

### Basic Usage

```python
from utils import read_tiff, calculate_tile_size
from calibration import estimate_calibration_from_meta, radiometric_calibration_lt1_with_linear
from rpc import RPCModel
from geometric_correction import ortho_rectify_by_rpc

# 1. Read SAR image and metadata
image_data, profile = read_tiff("input.tiff")
cal_params, product_level = estimate_calibration_from_meta("metadata.xml")

# 2. Radiometric calibration
sigma0_linear, sigma0_db = radiometric_calibration_lt1_with_linear(
    intensity_data,
    cal_params,
    polarization="VV"
)

# 3. Geometric correction with RPC
rpc_model = RPCModel.from_xml("rpc.xml")
ortho_data, ortho_meta = ortho_rectify_by_rpc(
    "sigma0.tif",
    rpc_model,
    output_path="ortho_sigma0.tif",
    output_crs="EPSG:4326",
    pixel_size=10.0  # 10m resolution
)
```

### Complete Processing Pipeline

See `main_example.py` for a complete processing workflow:

```bash
python main_example.py
```

This script demonstrates:
1. Metadata extraction
2. Radiometric calibration (amplitude → intensity → σ₀)
3. RPC-based orthorectification
4. Output generation in GeoTIFF format

## Module Reference

### calibration.py

Radiometric calibration functions for converting raw SAR intensities to calibrated backscatter.

**Key Functions:**
- `radiometric_calibration_lt1()` - LT-1 calibration (dB output)
- `radiometric_calibration_lt1_with_linear()` - Returns both linear and dB values
- `estimate_calibration_from_meta()` - Extract calibration parameters from metadata
- `get_polarization_type()` - Parse polarization information

**Calibration Formula (LT-1):**
```
σ₀ (dB) = 10 × log₁₀(I × (q/32767)²) - K

where:
  I = intensity (DN²)
  q = quantization level
  K = calibration constant
```

### geometric_correction.py

RPC-based geometric correction and orthorectification.

**Key Functions:**
- `ortho_rectify_by_rpc()` - Main orthorectification function
- `_make_coord_grid_for_block()` - Generate coordinate grids
- `_transform_coordinates_batch()` - Efficient batch coordinate transformation

**Features:**
- DEM-based height correction (optional)
- Customizable output projection (EPSG codes)
- Memory-efficient tiled processing
- Supports both geographic (lat/lon) and projected coordinate systems

### rpc.py

RPC (Rational Polynomial Coefficients) model implementation.

**Key Classes:**
- `RPCModel` - RPC coefficient container and transformation functions

**Supported Formats:**
- XML metadata files
- Plain text RPC files

### meta.py

Metadata parsing and bilinear geographic mapping model.

**Key Classes:**
- `GeoMappingModel` - Bilinear transformation for approximate georeferencing

**Use Case:**
- Fallback when RPC files are not available
- Quick approximate geolocation

### pol_maps.py

Polarimetric processing utilities.

**Key Functions:**
- `amplitude_block()` - Compute amplitude from complex SAR data
- `intensity_block()` - Compute intensity from amplitude

### incidence.py

Local incidence angle calculation for terrain correction.

### utils.py

Utility functions for I/O and data management.

**Key Functions:**
- `read_tiff()` - Read GeoTIFF files
- `iter_tiles()` - Generate tile windows for block processing
- `calculate_tile_size()` - Determine optimal tile size based on image dimensions

## Processing Workflow

```
┌─────────────────┐
│  Raw SAR Data   │
│   (Complex/     │
│   Intensity)    │
└────────┬────────┘
         │
         ▼
┌─────────────────┐
│   Amplitude     │  amplitude_block()
│   Extraction    │
└────────┬────────┘
         │
         ▼
┌─────────────────┐
│   Intensity     │  intensity_block()
│   Calculation   │
└────────┬────────┘
         │
         ▼
┌─────────────────┐
│   Radiometric   │  radiometric_calibration_lt1()
│   Calibration   │  → σ₀ (linear & dB)
└────────┬────────┘
         │
         ▼
┌─────────────────┐
│   Geometric     │  ortho_rectify_by_rpc()
│   Correction    │  → Georeferenced GeoTIFF
└────────┬────────┘
         │
         ▼
┌─────────────────┐
│  Final Product  │
│   (σ₀ GeoTIFF)  │
└─────────────────┘
```

## Output Products

The processing pipeline generates the following products:

| Product | Description | Format |
|---------|-------------|--------|
| `01_amplitude.tif` | Amplitude image (√I) | GeoTIFF (float32) |
| `01_intensity.tif` | Intensity image (I) | GeoTIFF (float32) |
| `02_sigma0_linear.tif` | Calibrated σ₀ (linear) | GeoTIFF (float32) |
| `02_sigma0dB.tif` | Calibrated σ₀ (dB) | GeoTIFF (float32) |
| `03_ortho_rpc_sigma0_linear.tif` | Orthorectified σ₀ | GeoTIFF (float32, geocoded) |

## Configuration

Edit `main_example.py` to configure processing parameters:

```python
# Input files
ROOT = r"E:/LT1/"
SRC_TIFF = ROOT + "1.tiff"      # SAR image
META_XML = ROOT + "1_meta.xml"  # Metadata
RPC_XML = ROOT + "1.rpc"        # RPC file
DEM_TIF = None                  # Optional DEM

# Output settings
OUT_DIR = ROOT + "output/"

# Geometric correction parameters
output_crs = "EPSG:4326"        # WGS84 geographic
pixel_size = 10.0               # 10 meter resolution
```

## Performance Considerations

- **Memory Usage**: Tiled processing is used to handle large images efficiently
- **Tile Size**: Automatically calculated based on image dimensions (typically 512-2048 pixels)
- **Processing Time**: Varies with image size and tile size
  - Example: 10000×10000 image typically processes in 5-15 minutes
- **Disk I/O**: LZW compression is applied to output GeoTIFFs to reduce file size

## Coordinate Systems

### Supported Projections

Any EPSG-defined coordinate system, commonly:
- **EPSG:4326** - WGS84 Geographic (lat/lon)
- **EPSG:3857** - Web Mercator
- **EPSG:32601-32660** - WGS84 UTM Northern Hemisphere
- **EPSG:32701-32760** - WGS84 UTM Southern Hemisphere

### RPC Model Requirements

- RPC coefficients must be provided in XML or text format
- Coordinates are internally handled in WGS84 (lat/lon)
- Output can be reprojected to any specified CRS

## Troubleshooting

### Common Issues

**Issue**: `ImportError: GDAL module not found`
```bash
# Solution: Install GDAL with proper bindings
conda install gdal
# or
pip install GDAL==$(gdal-config --version)
```

**Issue**: `Memory error during processing`
```python
# Solution: Reduce tile size in utils.calculate_tile_size()
TILE_SIZE = 512  # Reduce from default
```

**Issue**: `RPC file not found or invalid`
```python
# Solution: Verify RPC file format
rpc = RPCModel.from_xml("rpc.xml")  # For XML
# or
rpc = RPCModel.from_group_text_file("rpc.rpc")  # For text
```

**Issue**: `Geometric correction produces shifted results`
- Verify RPC coefficients are correct
- Check that DEM (if used) is in correct coordinate system
- Ensure pixel size is appropriate for the sensor resolution

## Data Formats

### Input Requirements

- **SAR Image**: GeoTIFF or similar raster format readable by `rasterio`
- **Metadata**: XML file containing calibration parameters and sensor info
- **RPC File**: XML or text file with RPC coefficients (optional for Level-1 products)
- **DEM** (optional): GeoTIFF in any projection (auto-reprojected internally)

### Metadata Structure (LT-1)

Expected XML structure:
```xml
<product>
  <generalHeader>
    <itemName>LT1A_...</itemName>
  </generalHeader>
  <acquisitionInfo>
    <polarisationMode>Single</polarisationMode>
  </acquisitionInfo>
  <imageDataInfo>
    <imageDataFormat>
      <calibrationConst K="..." q="..." />
    </imageDataFormat>
  </imageDataInfo>
</product>
```

## Contributing

Contributions are welcome! Please follow these guidelines:

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/new-feature`)
3. Commit your changes with clear messages
4. Add tests for new functionality
5. Submit a pull request

### Development Setup

```bash
# Clone repository
git clone https://github.com/yourusername/SARpreprocessing.git
cd SARpreprocessing

# Install development dependencies
pip install -r requirements.txt

# Run example
python main_example.py
```

## Citation

If you use this toolkit in your research, please cite:

```bibtex
@software{sar_preprocessing_toolkit,
  title = {SAR Preprocessing Toolkit},
  author = {Wu, Wenhao},
  year = {2025},
  url = {https://github.com/yourusername/SARpreprocessing}
}
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

- LT-1 satellite data format specifications
- RPC geometric correction algorithms from GDAL
- SAR radiometric calibration methods from ESA and CSA documentation

## Contact

- **Author**: Wu Wenhao
- **QQ**: 460249274
- **Issues**: Please report bugs and feature requests via [GitHub Issues](https://github.com/yourusername/SARpreprocessing/issues)

## References

### SAR Processing Literature

1. Cumming, I. G., & Wong, F. H. (2005). *Digital Processing of Synthetic Aperture Radar Data*. Artech House.
2. Moreira, A., et al. (2013). "A tutorial on synthetic aperture radar." *IEEE Geoscience and Remote Sensing Magazine*, 1(1), 6-43.
3. ESA Sentinel-1 Toolbox - Radiometric Calibration Documentation

### Related Tools

- **SNAP** (Sentinel Application Platform): ESA's official SAR processing software
- **GAMMA**: Commercial SAR/InSAR processing suite
- **ISCE** (InSAR Scientific Computing Environment): NASA/JPL InSAR processing framework
- **PyAR**: Python-based SAR processing library

---

**Version**: 1.0.0
**Last Updated**: December 2025
**Tested On**: Ubuntu 16.04.7 LTS, Python 3.7+
