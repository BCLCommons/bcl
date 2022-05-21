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

#ifndef BCL_MATH_FWD_HH_
#define BCL_MATH_FWD_HH_

// include bcl_defines.h header
#include "bcl_defines.h"

// external includes - sorted alphabetically
#include <functional>

// This file contains forward declarations for the math namespace
// This file is mostly automatically generated
// Developers may add typedefs and template default parameters
// all other changes will be removed the next time this file is generated
namespace bcl
{
  namespace math
  {
  /////////////////////
  // regular classes //
  /////////////////////

    class  Angle;
    class  AssignmentUnaryInterface;
    class  AssignmentUnaryStandard;
    class  BicubicSpline;
    class  ContingencyMatrix;
    class  ContingencyMatrixMeasures;
    class  CubicSpline;
    class  CubicSplineDamped;
    class  CubicSplineVariableDelta;
    class  DiscreteSetSelector;
    class  GaussianFunction;
    class  Gnuplot;
    class  GnuplotHeatmap;
    class  GnuplotMultiplot;
    class  Histogram;
    class  Histogram2D;
    class  Histogram3D;
    class  KernelFunction;
    struct LessThanAbsoluteMean;
    struct LessThanCounts;
    class  LinearFunction;
    class  LinearLeastSquares;
    class  LogLikelihood;
    class  MutateTransformationMatrix3D;
    class  MutateVector;
    class  PiecewiseFunction;
    class  Polynomial;
    class  QuadraticFunction;
    class  Quaternion;
    class  ROCCurve;
    class  RangeBorders;
    class  RotationMatrix2D;
    class  RotationMatrix3D;
    class  SmoothData;
    class  Statistics;
    class  TransformationMatrix3D;
    class  TricubicSpline;
    class  TrigonometricTransition;
    class  ZScore;

  //////////////////////
  // template classes //
  //////////////////////

    template< typename t_ArgumentType>
    class Assign;

    template< typename t_ArgumentType, template< typename> class t_ComparisonType>
    class AssignmentByComparison;

    template< typename t_ArgumentType>
    class AssignmentOperationInterface;

    template< typename t_DataType>
    class Assignments;

    template< typename t_ArgumentType1, typename t_ArgumentType2, typename t_ResultType>
    class BinaryFunctionBindFirst;

    template< typename t_ArgumentType1, typename t_ArgumentType2, typename t_ResultType>
    class BinaryFunctionBindSecond;

    template< typename t_ArgumentType1, typename t_ArgumentType2, typename t_ResultType>
    class BinaryFunctionCached;

    template< typename t_ArgumentType1, typename t_ArgumentType2, typename t_ResultType>
    class BinaryFunctionInterface;

    template< typename t_ArgumentType1, typename t_ArgumentType2, typename t_ResultType>
    class BinaryFunctionInterfaceSerializable;

    template< typename t_ArgumentType1, typename t_ArgumentType2, typename t_ResultType, typename t_Scalar>
    class BinarySumFunction;

    template< typename t_ArgumentType1, typename t_ArgumentType2, typename t_IntermediateType, typename t_ResultType>
    class BinaryUnaryFunctionAdapter;

    template< typename t_KeyType, typename t_KeyCompare = std::less< t_KeyType> >
    class Combination;

    template< typename t_DataType>
    class Comparisons;

    template< typename t_ArgumentType, typename t_ResultType>
    class ConstFunction;

    template< typename t_ArgumentType>
    class DivideEquals;

    template< typename t_ArgumentType, typename t_IntermediateType, typename t_ResultType>
    class FunctionAdapter;

    template< typename t_ArgumentType, typename t_ResultType>
    class FunctionCached;

    template< typename t_ArgumentType, typename t_ResultType>
    class FunctionInterface;

    template< typename t_ArgumentType, typename t_ResultType>
    class FunctionInterfaceSerializable;

    template< typename t_DataType>
    class Identity;

    template< typename t_ArgumentType>
    class MinusEquals;

    template< typename t_ArgumentType>
    class ModEquals;

    template< typename t_ArgumentType>
    class MutateCombine;

    template< typename t_DataType>
    class MutateDecisionNode;

    template< typename t_ArgumentType>
    class MutateInterface;

    template< typename t_ArgumentType>
    class MutateMoveWrapper;

    template< typename t_ArgumentType>
    class MutatePerturbation;

    template< typename t_ArgumentType>
    class MutateRepeat;

    template< typename t_ArgumentType>
    class MutateResult;

    template< typename t_ArgumentType, typename t_ResultType>
    class MutateTerminateDependent;

    template< typename t_ObjectType>
    class ObjectProbabilityDistribution;

    template< typename t_DataType>
    class ObjectStochasticSelector;

    template< typename t_ArgumentType>
    class PlusEquals;

    template< typename t_ArgumentType>
    class PowerEquals;

    template< typename t_DataType>
    class Range;

    template< typename t_DataType>
    class RangeSet;

    template< typename t_AverageType>
    class RunningAverage;

    template< typename t_AverageType>
    class RunningAverageSD;

    template< typename t_AverageType>
    class RunningMinMax;

    template< typename t_DataType>
    class RunningStatistics;

    template< typename t_SumType>
    class RunningSum;

    template< typename t_ArgumentType, typename t_ResultType>
    class SumFunction;

    template< typename t_Interface>
    class SumFunctionInterface;

    template< typename t_Interface>
    class SumFunctionMixin;

    template< typename t_DataType>
    class Tensor;

    template< typename t_ArgumentType>
    class TimesEquals;

  //////////////
  // typedefs //
  //////////////

  } // namespace math
} // namespace bcl

#endif // BCL_MATH_FWD_HH_
