# Mini Stack - SAR Interferometry Processing System

A comprehensive Synthetic Aperture Radar (SAR) interferometry processing framework designed for time-series InSAR analysis using the mini-stack approach. This system integrates with DORIS and StaMPS to perform Differential SAR Interferometry (DInSAR) and Persistent Scatterer (PS) analysis.

## Overview

Mini Stack is a Python-based processing pipeline that handles multi-temporal SAR data to detect ground deformation patterns. It implements an efficient "mini-stack" approach for processing large numbers of SAR image pairs by grouping them temporally and processing them in stages.

### Key Features

- **Mini-Stack Processing**: Efficient temporal grouping and processing of SAR interferograms
- **Differential InSAR**: Advanced DInSAR compression and analysis
- **Phase Linking**: Sophisticated phase linking techniques for improved accuracy
- **Statistical Analysis**: Kolmogorov-Smirnov testing for data quality assessment
- **DORIS Integration**: Seamless integration with DORIS processing workflows
- **StaMPS Support**: Full support for StaMPS PS analysis framework
- **Batch Processing**: Multi-threaded batch processing capabilities
- **Quality Filtering**: Advanced filtering mechanisms for reliable results

## Technology Stack

- **Language**: Python 2.x
- **Core Libraries**:
  - GDAL/OGR (geospatial data handling)
  - NumPy & SciPy (numerical computing)
  - scipy.optimize, scipy.stats, scipy.linalg
- **External Tools**:
  - DORIS (Delft object-oriented Radar Interferometric Software)
  - StaMPS (Stanford Method for Persistent Scatterers)
  - MATLAB (for StaMPS analysis)

## System Requirements

### Hardware Requirements
- **RAM**: Minimum 16GB (32GB+ recommended for large datasets)
- **Storage**: Adequate space for SAR data (typically 100GB+ per project)
- **CPU**: Multi-core processor recommended for parallel processing

### Software Requirements
- Python 2.7
- GDAL >= 2.0
- NumPy >= 1.10
- SciPy >= 0.17
- DORIS (latest stable version)
- StaMPS (version 4.x or higher)
- MATLAB (for StaMPS processing)

## Installation

### 1. Clone the Repository

```bash
git clone <repository-url>
cd mini_stack
```

### 2. Install Python Dependencies

```bash
# Install system dependencies (Ubuntu/Debian)
sudo apt-get update
sudo apt-get install python2.7 python-pip
sudo apt-get install gdal-bin libgdal-dev python-gdal

# Install Python packages
pip install numpy scipy
```

### 3. Install DORIS

```bash
# Download and install DORIS from the official repository
# Follow instructions at: https://github.com/TUDelftGeodesy/Doris
```

### 4. Install StaMPS

```bash
# Download StaMPS from: https://homepages.see.leeds.ac.uk/~earahoo/stamps/
# Follow the installation guide for your operating system
# Set environment variables:
export STAMPS=/path/to/stamps
export PATH=$PATH:$STAMPS/bin
```

### 5. Configure MATLAB (if using StaMPS MATLAB scripts)

Ensure MATLAB is installed and the StaMPS MATLAB scripts are in your MATLAB path.

## Project Structure

```
mini_stack/
├── Mini stack/              # Main processing pipeline
│   ├── mini_stack_prepare.py                    # Data preparation & scheduling
│   ├── insar2compress_input_crop_path.py        # Input processing
│   ├── mini_stack_compress_path_DSInSAR.py      # DInSAR compression
│   ├── mini_stack_time_compress_step_*.py       # Time-series compression stages
│   ├── concatenate_mini_stack_compress_minus_sign.py  # Concatenation
│   ├── Kol_smi_test.py                          # Statistical testing
│   └── dem.dorisin                              # DEM configuration
│
└── DS_Stamps/               # DORIS/StaMPS integration
    ├── DS_make.py                               # Main processing orchestrator
    ├── phase_linking_patch.py                   # Phase linking
    ├── shp_filter_patch.py                      # Quality filtering
    ├── DSInSAR_master_slave_combination.py      # Image combination
    ├── TOPS_Data_multilook.py                   # Multi-look processing
    └── ps_batch.m                               # MATLAB batch script
```

## Configuration

### 1. DEM Configuration

Edit the `dem.dorisin` file to configure your Digital Elevation Model:

```bash
cd "Mini stack"  # or DS_Stamps
nano dem.dorisin
```

Key parameters to configure:
- DEM file path
- DEM format (e.g., SRTM, ASTER)
- Memory allocation
- Processing verbosity

### 2. Input Data Preparation

Prepare your input data files:
- `Int_Data_merge.list` - List of interferogram pairs
- `day.1.in` or `slcs.list` - SLC image schedule
- Ensure SAR data is in MFF format with .res metadata


## run and build
FROM ubuntu:16.04

ENV DEBIAN_FRONTEND=noninteractive

# Install dependencies
RUN apt-get update && apt-get install -y \
    python2.7 \
    python-pip \
    gdal-bin \
    libgdal-dev \
    python-gdal \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

# Install Python packages
RUN pip install numpy scipy

# Copy project files
COPY . /app/mini_stack
WORKDIR /app/mini_stack

# Set environment variables
ENV PYTHONPATH=/app/mini_stack

CMD ["/bin/bash"]
```



- **Parallel Processing**: Utilize multi-core processors by adjusting batch processing parameters
- **Memory Management**: Monitor RAM usage and adjust processing chunk sizes
- **Storage**: Use SSD for faster I/O operations
- **Network**: For distributed processing, ensure high-bandwidth network connectivity



## Support

For issues and questions:
- Check the troubleshooting section
- Review DORIS documentation: https://github.com/TUDelftGeodesy/Doris
- Review StaMPS documentation: https://homepages.see.leeds.ac.uk/~earahoo/stamps/

## License

Please refer to the license file or contact the project maintainers for licensing information.

## Citation

If you use this software in your research, please cite:
- DORIS: https://github.com/TUDelftGeodesy/Doris
- StaMPS: Hooper, A., et al. (2012). "Recent advances in SAR interferometry time series analysis for measuring crustal deformation."
