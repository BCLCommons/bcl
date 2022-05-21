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

#ifndef BCL_ASSEMBLE_SSE_GEOMETRY_PACKER_INTERFACE_H_
#define BCL_ASSEMBLE_SSE_GEOMETRY_PACKER_INTERFACE_H_

// include the namespace header
#include "bcl_assemble.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_assemble_sse_geometry_interface.h"
#include "bcl_assemble_sse_geometry_packer_cache_wrapper.h"
#include "fold/bcl_fold_default_flags.h"
#include "math/bcl_math_binary_function_interface_serializable.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SSEGeometryPackerInterface
    //! @brief template class is the interface class for all SSEGeometry packers and pickers
    //! @details This interface provides access to minimal interface length and defines the operator that works on two
    //! SSEGeometryInterface derived classes. The return type is templated. It also provides the static function to
    //! wrap any SSEGeometryPackerInterface derived class with a cache function.
    //!
    //! @tparam t_ReturnType whether this packer returns a single Packing, a list or a list of lists or etc. of SSEGeometryPackings
    //!
    //! @remarks example unnecessary
    //! @author karakam
    //! @date Aug 24, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_ReturnType>
    class SSEGeometryPackerInterface :
      public math::BinaryFunctionInterfaceSerializable< SSEGeometryInterface, SSEGeometryInterface, t_ReturnType>
    {

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief Clone function
      //! @return pointer to new SSEGeometryPackerInterface< t_ReturnType>
      virtual SSEGeometryPackerInterface< t_ReturnType> *Clone() const = 0;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      //! @brief get minimal interface length
      //! @return get minimal interface length
      virtual double GetMinimalInterfaceLength() const = 0;

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief allows a SSEGeometryPackerInterface derived class to be wrapped in a cache function if the flag is provided
      //! @param SP_PACKER_FUNCTION ShPtr to SSEGeometryPackerInterface derived class of interest
      //! @param CACHEABLE whether this function should be cached if the flag is provided
      //! @param IS_SYMMETRIC whether the provided function is symmetric
      //! @return ShPtr to wrapped class
      static util::ShPtr< SSEGeometryPackerInterface< t_ReturnType> > WrapPacker
      (
        const util::ShPtr< SSEGeometryPackerInterface< t_ReturnType> > &SP_PACKER_FUNCTION,
        const bool CACHEABLE,
        const bool IS_SYMMETRIC
      )
      {
        // check if the packing list function is going to be cached
        if( CACHEABLE)
        {
          // warp function into cache function
          util::ShPtr< SSEGeometryPackerInterface< t_ReturnType> > sp_cache_function
          (
            new SSEGeometryPackerCacheWrapper< t_ReturnType>( SP_PACKER_FUNCTION, IS_SYMMETRIC)
          );

          // insert and return the inserted
          return sp_cache_function;
        }
        // if not caching is necessary
        else
        {
          return SP_PACKER_FUNCTION;
        }
      }

    }; // template class SSEGeometryPackerInterface

  } // namespace assemble
} // namespace bcl

#endif // BCL_ASSEMBLE_SSE_GEOMETRY_PACKER_INTERFACE_H_
