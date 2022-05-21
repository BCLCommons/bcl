// (c) Copyright BCL @ Vanderbilt University 2014
// (c) BCL Homepage: http://www.meilerlab.org/bclcommons
// (c) BCL Code Repository: https://github.com/BCLCommons/bcl
// (c)
// (c) The BioChemical Library (BCL) was originally developed by contributing members of the Meiler Lab @ Vanderbilt University.
// (c)
// (c) The BCL is now made available as an open-source software package distributed under the permissive MIT license,
// (c) developed and maintained by the Meiler Lab at Vanderbilt University and contributing members of the BCL Commons. 
// (c)
// (c) External code contributions to the BCL are welcome. Please visit the BCL Commons GitHub page for information on how you can contribute.
// (c) 
// (c) This file is part of the BCL software suite and is made available under the MIT license.
// (c)

#ifndef BCL_OPENCL_FWD_HH_
#define BCL_OPENCL_FWD_HH_

// include the dependency file for this header
#include "bcl_opencl.depends.fwd.hh"

// This file contains forward declarations for the opencl namespace
// This file is mostly automatically generated
// Developers may add typedefs and template default parameters
// all other changes will be removed the next time this file is generated
namespace bcl
{
  namespace opencl
  {
  /////////////////////
  // regular classes //
  /////////////////////

    class ApproximatorResilientPropagation;
    class ApproximatorSequentialMinimialOptimization;
    class ApproximatorSimplePropagation;
    class Buffer;
    class CommandQueue;
    class Context;
    class DensityFitProteinMinimzerPowell;
    class DensitySimulateGaussianSphere;
    class Device;
    class ExtensionData;
    class Extensions;
    class InsertionSort;
    class KappaNearestNeighbor;
    class KernelSourceAlternative;
    class KernelSourceFile;
    class KernelSourceInterface;
    class KernelSourceString;
    class KernelSources;
    class ModelInterface;
    class NeuralNetwork;
    class Platform;
    class ProteinAgreementCCC;
    class QualityGDT;
    class QualityLCS;
    class QualityRMSD;
    class RMSD;
    class SaxsDebye;
    class SupportVectorMachine;
    class Tools;

  //////////////////////
  // template classes //
  //////////////////////

    template< typename t_DataType>
    class ArgMax;

    template< typename t_DataType>
    class ArgMin;

    template< typename t_DataType>
    class CoordinateTransformer;

    template< typename t_DataType>
    class DataSetMinMax;

    template< typename t_DataType>
    class EuclideanDistance;

    template< typename t_DataType>
    class FeatureSimilarityMeasures;

    template< typename t_DataType>
    class Matrix;

    template< typename t_DataType>
    class Matrix3x3;

    template< typename t_DataType>
    class MatrixAdd;

    template< typename t_DataType>
    class MatrixMultiply;

    template< typename t_DataType>
    class MatrixTranspose;

    template< typename t_DataType>
    class Operations;

    template< typename t_DataType>
    class SingularValueDecomposition;

    template< typename t_DataType>
    class TransferFunctionGaussian;

    template< typename t_DataType>
    class TransferFunctionSigmoid;

    template< typename t_DataType>
    class Vector;

    template< typename t_DataType>
    class VectorMatrixAdd;

  //////////////
  // typedefs //
  //////////////

    typedef util::Enum< ExtensionData, Extensions>                          Extension;
    typedef util::Enum< util::ShPtr< KernelSourceInterface>, KernelSources> KernelSource;

  } // namespace opencl
} // namespace bcl

#endif // BCL_OPENCL_FWD_HH_
