<a href="http://www.physycom.unibo.it"> 
<div class="image">
<img src="https://cdn.rawgit.com/physycom/templates/697b327d/logo_unibo.png" width="90" height="90" alt="© Physics of Complex Systems Laboratory - Physics and Astronomy Department - University of Bologna"> 
</div>
</a>
<a href="https://travis-ci.org/physycom/inertial_tools"> 
<div class="image">
<img src="https://travis-ci.org/physycom/inertial_tools.svg?branch=master" width="90" height="20" alt="Build Status"> 
</div>
</a>
<a href="https://ci.appveyor.com/project/cenit/inertial-tools"> 
<div class="image">
<img src="https://ci.appveyor.com/api/projects/status/c6vaxn3m79r82y2m?svg=true" width="90" height="20" alt="Build Status"> 
</div>
</a>

## Purpose
This document describes a collection of tools, developed within a standardized framework, to perform a number of quantitative and qualitative analysis of inertial data generated by various sources. With the term inertial we refer mainly to those dataset which describes the state of motion of an ideal rigid body, they are:
- a 3-dimensional vector containing the rigid body linear acceleration;
- a 3-dimensional vector containing the rigid body angular velocity;
- a scalar value representing the modulus of the linear velocity vector;
- an instant of time to which all the previous quantities are refered, denoted as _timestamp_.

Starting from the knowledge of a time series of inertial data it is then possible to perform several analysis which ranges from 2D/3D motion reconstruction for small time intervals (8-10 seconds) to the estimation of statistical quantities which provides a number of markers describing the motion itself and the global quality of the dataset. The statistical analysis are usually performed on large time scales, from hours up to weeks.

The suite is composed of tools which match a specific analysis pattern/method and are intended to serve a specific purpose, as described below. The tools are:
- `csv_to_inertial`:

	As we discovered by facing dataset coming from different sources (vendors, models, firmwares, devices) the development of a standardized format is a milestone. This tool allows to convert from a generic CSV dataset to our `inertial` [format](https://github.com/physycom/file_format_specifications). 

- `data_calibrator`:

	This tool is a workstation implementation of the _Continuous Autonomous Calibration Algorithm_, see the documentation of [libcasc](https://github.com/physycom/libcasc) for more details. The precise purpose of this tool is to constitute a benchmark utility which has been used to fine-tune the algorithm's parameters along with its performances.

- `data_filter`

	Due to the usually high rate of output of common inertial sensors (from 2 kHz down to 200 Hz) the amount of data produced by these types of devices is very large. Since most of algorithms are usually concerned by a specific combination of ranges of values for the physical quantities. This tool provides a simple and fast utility to extract a subset of a given dataset, tipically large, to reduce the dimensions and increase performances. The aim of such an operation is valuable when one is faced with the problem of optimizing and fine-tuning various purpose algorithms.

- `data_generator`

	This is a support tool which somehow reverts the process of rigig body motion reconstruction. Starting from a given (hardcoded) spatial trajectory for the rigid body's center of mass this tool creates a discretized version of the inertial quantities (in the sense precised above) associated to such a motion. This tool has revealed crucial to test advanced motion reconstrution algorithm since the possibility of creating a dataset from a fully analytically controlled motion offers the opportunity to test, fine-tune and debug in a supervised way.

	*IMMAGINE DEI DATI GENERATI PER LA CHICANE*

- `data_purge`

	As a matter of fact the data filtering operation is crucial when one is addressed to numerical dataset analysis. In our case, due to the intrinsic noisyness of the inertial data coming from MEMS cheap sensors, the filtering operation is even more difficult since the information carried by our data spans a very broad frequency band and are therefore difficult to detect. The quest for optimal filter parameter is made even harder by the trial-and-error intrinsic nature of the process. To overcome this difficulty we developed a flexible and performant filtering program which is capable of applying a user-defined filter pattern to a given dataset. With the aid of this tool it is then possible to automatize to process of optimizing the filter shape and parameters, increasing the effectiveness of the filter itself.

	*IMMAGINE DELLA PURGA*

- `data_reconstruction`

	Motion reconstruction from MEMS inertial data is a challenging task. The extreme noisiness of the datasets is difficult to control numerically and the extent of the reconstruction is confined within several seconds (8-10) from the initial condition, which is as well difficult to synchronize due to the fact that is usually provided by a different sensor. Nonetheless we managed to collect various previous tool developed by group into a single standardized tool which perform a numerical integration of the inertial data to reconstruct the motion of the rigid body.

	*QUALCOSA DI RICOSTRUITO*

- `data_rotator`

	A tool which enables the manipulation of a given dataset in terms of rotation. It has been developed as an intermediate layer between the original dataset and the `data_reconstruction` tool with aim of referring the inertial data to the correct frame of reference, i.e x axis in the direction of motion, y axis orthogonal to it and a vertical z axis.

	*IMMAGINE DI DATI RUOTATI*

- `inertial_stats`

	A numerical tool to analyze a given dataset in terms of the estimation of statistical quantities, such as:
	- mean values, for the inertial data;
	- quadratic moments and standard deviation;
	- quadratic covariances and cross-correlations.

	These numerical quantities are useful for a number of purposes: from a quantitative noise estimation of the device to the dynamical statistical characterization of a given dataset, in terms of driving behavior or overall device functionality monitoring.

	*IMMAGINE DI NOISE TEST*

- `mutual_entropy`

	This tool implements the so called Mutual Information technique to estimate numerically the degree of coeherence between the different inertial quantities. By employing the evaluation of a parameter called _Mutual Entropy_ it is possible to classify quantitatively the performance of accelerometer and gyroscope related to the same dataset.

	*IMMAGINE MUTUAL ENTROPY*

- `speed_compare`

	This tool provides a quantitative estimate of the time coherence and synchronicity between the inertial sensors and the GNSS module. This is achieved by implementing an integration scheme which is capable of reconstructing the modulus of the linear velocity vector in the two-dimensional plane of motion of the device, considering only acceleration data. Once this is done a smoothing and filtering process on the velocity data is performed (we recall that in usual situation the rate of the inertial data is roughly an order of magnitude greater than the rate of the geopositioning data) and then the two velocity dataset are compared with respect to time. The purpose of this tool is twofold:
	- estimate the given degree of coherence and synchronicity;
	- provide a quantitative way to correct, in post-processing, any discrepancy with the aim of aiding the `data_reconstruction` algorithm.

	*IMMAGINE SPEED COMPARE, QUESTA BUONA MAI AVUTA*

This suite has been developed by the Physics of the City laboratory, Physics and Astronomy Department, University of Bologna, as an internal support tool to support the consultancy activity within the official collaboration with UnipolSai.

## Technical notes
The repository comes equipped with both a *makefile* and a *VS15* solution. 

Every tool of the suite provides its own _synopsis_ when run on a terminal.

Most of the tools produce a graphical result, constituted by a `gnuplot` script, thus to enjoy this feature gnuplot is required.
