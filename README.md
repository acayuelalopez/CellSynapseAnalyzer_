# Cell Synapse Analyzer Plugin for ImageJ

**Cell Synapse Analyzer** is a plugin for ImageJ designed to facilitate the analysis of cell synapses in microscopy images. It provides a graphical interface for selecting, processing, and analyzing synapse data, focusing on the quantification of synaptic connections and their spatial distribution.

## Main Features

- Load and preview multiple images from a directory.
- Select and process images to identify and quantify synaptic connections.
- Automatic and manual modes for synapse analysis.
- Calculation of synapse density and spatial distribution.
- Export results in CSV format.
- Save processed images as TIFF files.

## Requirements

- Java 8 or higher.
- ImageJ (preferably the Fiji distribution).
- Apache Commons Math library for statistical analysis.
- JFreeChart library for plotting.

## Installation

1. Compile or download the JAR file of the plugin.
2. Place it in the `plugins` folder of ImageJ/Fiji.
3. Restart ImageJ.
4. Access the plugin via `Plugins > Cell Synapse Analyzer`.

## Workflow

1. Select a directory of images.
2. Define a region of interest (ROI) with the rectangle tool.
3. Process the image to identify and quantify synaptic connections.
4. View and adjust the ROIs if necessary.
5. Export the results and save the processed images.

## Outputs

- CSV file with quantitative data per image.
- TIFF images of the masks and generated montages.
- Individual and total results tables.

## Application

This plugin is useful for researchers in neuroscience and cell biology, enabling semi-automated quantification of synaptic connections in microscopy images.

## License

Distributed under the MIT License.

