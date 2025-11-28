# PyPhotonicReagents

A standalone Pythonic suite for advanced laser experiments involving a broad class of equipment, such as laser pulse shapers, spectrometers, cameras, and moving stages.

## Overview
Note: This is implemented in Python 2.
PyPhotonicReagents provides a unified Python interface for controlling laboratory instruments commonly used in ultrafast laser experiments and photonics research. Developed by [Prof. Denys Bondar](https://sse.tulane.edu/denys-bondar-phd) at Tulane University, this library supports experimental work in quantum control and ultrafast nonlinear optics.

## Supported Equipment

- **Laser Pulse Shapers** - Spatial Light Modulator (SLM) based devices for phase and amplitude modulation of ultrafast pulses
- **Spectrometers** - Optical spectrum analyzers for spectral characterization
- **Cameras** - Scientific imaging devices for beam profiling and data acquisition
- **Moving Stages** - Motorized translation/rotation stages for optical alignment and scanning

## Installation

```bash
git clone https://github.com/dibondar/PyPhotonicReagents.git
cd PyPhotonicReagents
pip install -r requirements.txt
```

## Peripheral Interaction Architecture

PyPhotonicReagents implements peripheral communication through a modular wrapper-based architecture. This design pattern enables:

### Communication Layer Design

1. **Vendor SDK Wrappers**
   - Each hardware device communicates through manufacturer-provided SDKs (DLLs, shared libraries, or serial protocols)
   - Python wrappers encapsulate low-level calls using `ctypes` or vendor-specific Python bindings
   - This abstraction isolates device-specific quirks from experimental logic

2. **Hardware Abstraction Layer**
   - Uniform API across different device types (e.g., all spectrometers share common `acquire()`, `get_spectrum()` methods)
   - Device configuration handled through standardized parameter dictionaries
   - Enables swapping equipment without modifying experiment scripts

3. **Connection Protocols**
   - **USB/Serial**: For spectrometers and stage controllers via PySerial or direct USB
   - **HDMI/DVI Display**: Pulse shapers using SLMs typically appear as secondary monitors; phase masks are rendered as images
   - **Ethernet/TCP-IP**: Network-connected devices for remote control
   - **Camera Interfaces**: Vendor SDKs (e.g., DCAM, uEye) or standards like GenICam



### Synchronization and Triggering

Laboratory experiments often require precise timing between components:

- **Software Triggering**: Python controls acquisition sequences directly
- **Hardware Triggering**: TTL signals coordinate laser shots with data acquisition
- **Event-Driven Callbacks**: Some SDKs support asynchronous notification of acquisition completion


## Related Research

This library supports research in:
- Coherent quantum control of chemical reactions
- MIIPS (Multiphoton Intrapulse Interference Phase Scan) pulse compression
- Optimal Dynamic Discrimination (ODD) for fluorescent protein characterization
- Ultrafast spectroscopy and pulse characterization

## Author

**Denys I. Bondar, Ph.D.**  
- [ORCID](https://orcid.org/0000-0002-3626-4804)



> **Note**: This README was generated based on the repository description and author's research context. For implementation details, refer to the source code in the repository.
