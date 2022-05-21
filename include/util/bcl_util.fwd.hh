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

#ifndef BCL_UTIL_FWD_HH_
#define BCL_UTIL_FWD_HH_

// include the dependency file for this header
#include "bcl_util.depends.fwd.hh"

// This file contains forward declarations for the util namespace
// This file is mostly automatically generated
// Developers may add typedefs and template default parameters
// all other changes will be removed the next time this file is generated
namespace bcl
{
  namespace util
  {
  /////////////////////
  // regular classes //
  /////////////////////

    class Assert;
    class CPPDataTypes;
    class CallStack;
    class ClassNameStandardizer;
    class CleanableInterface;
    class ColorGradient;
    class Colors;
    class DataType;
    class EnumsInstances;
    class Format;
    class FunctionalType;
    class ImplementationInterface;
    class LoggerDefault;
    class LoggerInterface;
    class Loggers;
    class MemoryUsage;
    class Message;
    class ObjectDataLabel;
    class ObjectDataLabelTokenizer;
    class ObjectInstances;
    class ObjectInterface;
    class RuntimeEnvironmentDefault;
    class RuntimeEnvironmentInterface;
    class RuntimeEnvironments;
    class SerializableInterface;
    class Stopwatch;
    class StringReplacement;
    class Time;
    class UndefinedObject;

  //////////////////////
  // template classes //
  //////////////////////

    template< typename t_ArgumentType1, typename t_ArgumentType2, typename t_ResultType>
    class BinaryFunctionInterface;

    template< typename t_ArgumentType1, typename t_ArgumentType2, typename t_ResultType>
    class BinaryFunctionInterfaceNonConst;

    template< typename t_ArgumentType1, typename t_ArgumentType2, typename t_ResultType>
    class BinaryFunctionInterfaceSerializable;

    template< typename t_BinaryFunction>
    class BinaryFunctionSTLWrapper;

    template< typename t_ArgumentType1, typename t_ArgumentType2, typename t_ResultType, typename t_FunctionClass>
    class BinaryFunctionWrapper;

    template< typename t_DataType>
    class CPUBenchmarkWhetstone;

    template< typename t_DataType, typename t_Derived>
    class Enum;

    template< typename t_DataType>
    class EnumData;

    template< typename t_DataType, typename t_Derived>
    class Enumerate;

    template< typename t_Interface>
    class Enumerated;

    template< typename t_ArgumentType1, typename t_ArgumentType2, typename t_ArgumentType3, typename t_ArgumentType4, typename t_ResultType>
    class FourInputFunctionInterface;

    template< typename t_ArgumentType, typename t_ResultType>
    class FunctionInterface;

    template< typename t_ArgumentType, typename t_ResultType>
    class FunctionInterfaceNonConst;

    template< typename t_ArgumentType, typename t_ResultType>
    class FunctionInterfaceSerializable;

    template< typename t_ArgumentType, typename t_ResultType, typename t_FunctionClass>
    class FunctionWrapper;

    template< typename t_Interface>
    class Implementation;

    template< typename t_DataType>
    class OwnPtr;

    template< typename t_DataType>
    class PtrInterface;

    template< typename t_DataType>
    class ShPtr;

    template< typename t_DataType>
    class ShPtrList;

    template< typename t_DataType>
    class ShPtrVector;

    template< typename t_DataType>
    class SiPtr;

    template< typename t_DataType>
    class SiPtrList;

    template< typename t_DataType>
    class SiPtrVector;

    template< typename t_ArgumentType1, typename t_ArgumentType2, typename t_ArgumentType3, typename t_ResultType, typename t_FunctionClass>
    class TertiaryFunctionWrapper;

    template< typename t_ArgumentType1, typename t_ArgumentType2, typename t_ArgumentType3, typename t_ResultType>
    class ThreeInputFunctionInterface;

    template< typename t_ResultType>
    class ThunkInterface;

    template< typename t_FunctionClass, typename t_ResultType>
    class ThunkWrapper;

    template< typename t_DataType>
    class VoxelGrid;

    template< typename t_DataType>
    class Wrapper;

    template< typename t_BaseType>
    class WrapperBase;

    template
    <
      typename t_Enum,
      const std::string &( *GetEnumDescriptor)( const t_Enum &),
      size_t s_NumberOfEnums
    >
    class WrapperEnum;

  //////////////
  // typedefs //
  //////////////

    typedef Enum< linal::Vector3D, Colors>                                  Color;
    typedef Enum< ShPtr< LoggerInterface>, Loggers>                         Logger;
    typedef Enum< ShPtr< RuntimeEnvironmentInterface>, RuntimeEnvironments> RuntimeEnvironment;

  } // namespace util
} // namespace bcl

#endif // BCL_UTIL_FWD_HH_
