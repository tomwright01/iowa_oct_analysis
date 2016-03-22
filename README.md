# IOWA OCT Analysis

A python class for manipulating optical coherence tomography (OCT) segmentation data from the [Iowa Reference algorithm](https://www.iibi.uiowa.edu/content/iowa-reference-algorithms-human-and-murine-oct-retinal-layer-analysis-and-display).

The Iowa reference algorithm uses a graph based approach to identify 10 surfaces in retinal OCT's. The GUI interface enables visualisation of retinal layers and can calculate retinal thickness values between any pair of surfaces. Thickness values can be summarised within (amonst others) the 9 ETDRS retinal regions.

The GUI does not support batch processing of multiple OCT's or export of summarised thickness values in a machine readable format. This class tries to bridge that gap, enabling batch processing and automated analysis of data from multiple OCT recordings.

The class also supports loading raw OCT data enabling the surface information to be overlayed back onto the original OCT. This makes extracting layer intensity information and producing graphics easy.

## Segmenting OCT recordings

The Iowa reference algorithm GUI is used for the actual segmentation. Currently the surfaces and grid center data are supported.
Minimal export required:
* XXX_Surfaces_Iowa.xml
* XXX_GridCenter_Iowa.xml

## Raw OCT data

Currently only data from the Cirrus OCT system (Carl Zeiss Meditec AG) is supported. Raw OCT data can be exported in .img format. This requires the Research Bundle (available under separate licence) to be installed.
Support for other systems is planned.

## Example usage

See file `examples.py` for example use.

## Licence
This work is licensed under the terms of the MIT license.
