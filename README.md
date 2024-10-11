# OKM Research Project

## Overview

This project implements various k-means clustering algorithms in C++ to analyze and process different datasets. The objective is to explore the efficiency and accuracy of batch k-means, online k-means, and their incremental variants. By applying these algorithms to multiple datasets, the project evaluates their performance based on metrics such as Sum of Squared Errors (SSE) and computation time.

## Features

- **Data Management**: Structures and functions to handle datasets, including dynamic memory allocation, loading from files, and normalization.
- **K-Means Variants**: Implements Batch K-Means (BKM), Online K-Means (OKM), and Incremental variants like IBKM and IOKM.
- **Initialization Methods**: Supports K-Means++ initialization, random selection, and Maximin for initializing cluster centers.
- **Performance Measurement**: Calculates SSE and measures execution time for each algorithm variant.
- **Replication and Statistics**: Performs multiple runs to gather statistics on SSE and iteration counts, providing insights into algorithm stability and performance.
- **Modular Design**: Clear separation of data structures, memory management, and algorithm implementations for ease of maintenance and extension.

## Detailed Description

### Data Structures

- **Data_Set**: Represents a dataset with attributes including data points, number of points, dimensions, categories, and category labels.

- **Run**: Stores results from multiple runs of clustering algorithms, including initial and final SSE values and the number of iterations for each run.

- **Dbl_Stats and Int_Stats**: Structures to store statistical data such as minimum, maximum, mean, and standard deviation for double and integer metrics respectively.

### Memory Management

Functions are provided to allocate and free memory for datasets and runs. This includes allocating 2D arrays for data points and managing memory for statistical data.

### Data Handling

- **Loading Data**: The load_data_set function reads datasets from files, expecting a specific format where the first line contains the number of points, dimensions, and categories, followed by data points with optional category labels.

- **Normalization**: The min_max_normalize function normalizes dataset attributes to a [0,1] range using min-max normalization, which is crucial for ensuring that all attributes contribute equally to the distance calculations.

- **Resetting Data**: Functions like reset_data_set reset the dataset by setting all attribute values and categories to zero, preparing the data structure for new clustering assignments.

### Distance Calculation and SSE

- **Distance Function**: The sqr_euc_dist function calculates the squared Euclidean distance between two vectors, a fundamental operation in k-means clustering.

- **Sum of Squared Errors (SSE)**: The calc_sse function computes the SSE for a given clustering by summing the squared distances of each point to its nearest cluster center.

### K-Means Algorithms

- **Batch K-Means (BKM)**: This implementation assigns each data point to the nearest cluster center and updates the centers based on the mean of assigned points. The process repeats until no changes in assignments occur.

- **Online K-Means (OKM)**: Updates cluster centers incrementally as points are processed one by one, using a learning rate to adjust the centers based on new data.

- **Incremental Batch K-Means (IBKM)** and **Incremental Online K-Means (IOKM)**: These variants implement incremental approaches to k-means clustering by splitting clusters and refining centers in a stepwise manner.

### Initialization Methods

- **K-Means++**: The kmeanspp function implements the K-Means++ initialization method, which selects initial cluster centers with a probability proportional to the squared distance from existing centers to enhance clustering performance.

- **Random Selection**: The rand_sel function randomly selects initial cluster centers from the dataset, providing a baseline initialization method.

- **Maximin**: The maximin function selects initial cluster centers by choosing points that are maximally distant from existing centers, ensuring diverse initial centers.

### Replication and Statistics

- **Replication Functions**: Functions like rep_rand_sel_okm, rep_maximin_okm, and rep_kmeanspp_okm perform multiple runs of clustering algorithms with different initialization methods to gather statistical data.

- **Statistics Calculation**: The comp_stats function computes statistical measures such as minimum, maximum, mean, and standard deviation for SSE and iteration counts across multiple runs, offering insights into algorithm performance and variability.

### Main Function

The main function orchestrates the execution of clustering algorithms on a series of datasets:

1. **Dataset Selection**: Iterates over a predefined list of dataset files and corresponding cluster counts.

2. **Clustering Execution**:
   - Applies Batch K-Means with K-Means++ initialization and measures SSE and execution time.
   - Applies Online K-Means with K-Means++ initialization and measures SSE and execution time.
   - Applies Incremental Batch K-Means and Online K-Means followed by Batch K-Means, measuring performance metrics.

3. **Performance Reporting**: Outputs the final SSE and execution time for each clustering algorithm applied to each dataset.

4. **Memory Cleanup**: Frees allocated memory after processing each dataset to prevent memory leaks.

### Performance Measurement

The code utilizes the std::chrono library to measure the execution time of clustering algorithms in milliseconds, providing quantitative metrics for evaluating the efficiency of each method.

### Random Number Generation

The project employs the Mersenne Twister std::mt19937 for random number generation, ensuring reproducibility and high-quality randomness in initialization methods like K-Means++ and random selection.

### Error Handling

The code includes error checking for memory allocation and data loading, ensuring that the program exits gracefully with informative messages if issues arise, such as insufficient data or incorrect file formats.

## Conclusion

The OKM Research Project offers a robust and flexible implementation of various k-means clustering algorithms, enabling comprehensive performance analysis across diverse datasets. By incorporating multiple initialization strategies and algorithmic variations, the project facilitates a deep understanding of k-means clustering dynamics, efficiency, and accuracy. The modular design and detailed statistical analysis support further research and optimization in clustering methodologies.
