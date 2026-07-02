# TMS-EMG Kit

*Version 2.0*

The TMS-EMG ToolBox is a data analysis application built for researchers working with Transcranial Magnetic Stimulation (TMS) and electromyography (EMG) data. It takes raw, multi-channel recordings from a wide range of common data acquisition systems and converts them into one standardized format for analysis.

Because different labs often use different recording hardware, comparing or combining datasets can be a major hassle. This toolbox solves that problem: whether your data comes from a brand-new experiment or a retrospective study collected on entirely different equipment, it can be imported, organized, processed, and analyzed the same way. That makes it much easier for research sites to pool data and collaborate.

### Supported Import Formats

- Labchart (.mat)
- Brainsight (.txt)
- Signal (.mat, .cfs)
- Acqknowledge (.mat, .acq)
- BrainVision (.eeg)
- Spike (.smr)
- Custom file types

## Getting Started

The toolbox is available in two forms:

- **MATLAB-Independent** — a standalone application that does not require a MATLAB license.
- **MATLAB-Dependent** — runs inside MATLAB (Windows and Mac compatible).

**Requirements for the MATLAB-Dependent version:**
- MATLAB 2023a
- Statistics and Machine Learning Toolbox
- Signal Processing Toolbox
- Mapping Toolbox

### Quick Start

1. Download the latest release from this repository.
2. Install MATLAB 2023a and the required toolboxes listed above (skip this step if you're using the MATLAB-Independent version).
3. Launch the application.
4. On first launch, you'll be prompted to select a data folder (where your raw files live) and a save folder (where results are stored).
5. Follow the full tutorial for a walkthrough of importing, organizing, processing, and analyzing your data: https://cunninghamlab.framer.website/tools

## Examples

**Full Overview**

![Full Overview](Full%20Overview.jpg)

**Data Organization / Epoching**

![Epoch](Epoch.jpg)

**Signal Processing**

![Signal Processing](SignalProcessing.jpg)

**MEP Analysis**

![MEP Analysis](MEP%20Analysis.jpg)

**Silent Period Analysis**

![Silent Period Analysis](Silent%20Period%20Analysis.jpg)

**Export Outcomes and Settings (Long Format)**

![Export Outcomes and Settings Long Format](Export%20Outcomes%20and%20Settings%20Long%20Format.jpg)

## Stay Up to Date

Join our Google Group for announcements and updates: https://groups.google.com/a/case.edu/g/tmsemgkit
