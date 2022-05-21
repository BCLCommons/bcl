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

#ifndef BCL_IO_FWD_HH_
#define BCL_IO_FWD_HH_

// include the dependency file for this header
#include "bcl_io.depends.fwd.hh"

// This file contains forward declarations for the io namespace
// This file is mostly automatically generated
// Developers may add typedefs and template default parameters
// all other changes will be removed the next time this file is generated
namespace bcl
{
  namespace io
  {
  /////////////////////
  // regular classes //
  /////////////////////

    class BinarySerialize;
    class Directory;
    class DirectoryEntry;
    class File;
    class FileStreamBuffer;
    class FixedLineWidthWriter;
    class IFStream;
    class OFStream;
    class Serialization;
    class SerializationInterface;
    class Serialize;
    class Serializer;
    class StreamBufferClasses;
    class StreamBufferInterface;
    class StreamInterface;
    class ValidationResult;

  //////////////////////
  // template classes //
  //////////////////////

    template< typename t_Type, typename t_EnsembleType>
    class RetrieveInterface;

    template< typename t_DataType>
    class SerializationBase;

    template< typename t_DataType>
    class SerializationBuiltin;

    template< typename t_ContainerType>
    class SerializationContainer;

    template< typename t_ContainerType>
    class SerializationMap;

    template< typename t_DataType>
    class SerializationViaStaticFunctions;

    template< typename t_DataType>
    class SerializationWithCheck;

    template< typename t_DataType>
    class SerializationWithMinMax;

    template< typename t_DataType>
    class SerializationWrapper;

    template< typename t_StoredType>
    class StoreInterface;

  //////////////
  // typedefs //
  //////////////

    typedef util::Enum< util::ShPtr< StreamBufferInterface>, StreamBufferClasses> StreamBufferClass;

  } // namespace io
} // namespace bcl

#endif // BCL_IO_FWD_HH_
