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

#ifndef BCL_ASSEMBLE_SSE_GEOMETRY_PACKER_CACHE_WRAPPER_H_
#define BCL_ASSEMBLE_SSE_GEOMETRY_PACKER_CACHE_WRAPPER_H_

// include the namespace header
#include "bcl_assemble.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_assemble_sse_geometry_packer_interface.h"
#include "math/bcl_math_binary_function_cached.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SSEGeometryPackerCacheWrapper
    //! @brief This class allows wrapping a SSEGeometryPackerInterface derived class with a cache function
    //! @details This template class allows wrapping a packer class with cache while allowing to keep the wrapped function still as a
    //! SSEGeometryPackerInterface derived class and thus have access to minimal interface length used.
    //!
    //! @tparam t_ReturnType Return Type of the SSEGeometryPackerInterface derived class
    //!
    //! @see @link example_assemble_sse_geometry_packer_cache_wrapper.cpp @endlink
    //! @author karakam
    //! @date Aug 24, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_ReturnType>
    class SSEGeometryPackerCacheWrapper :
      public SSEGeometryPackerInterface< t_ReturnType>
    {

    private:

    //////////
    // data //
    //////////

      //! minimal interface length
      double m_MinimalInterfaceLength;

      //! cache function
      util::ShPtr
      <
        math::BinaryFunctionCached< SSEGeometryInterface, SSEGeometryInterface, t_ReturnType>
      > m_CacheFunction;

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      SSEGeometryPackerCacheWrapper< t_ReturnType>() :
        m_MinimalInterfaceLength(),
        m_CacheFunction()
      {
      }

      //! @brief constructor from a SSEPackerInterface derived class to wrap
      //! @param SP_PACKER_FUNCTION ShPtr to Packer function to wrap
      //! @param IS_SYMMETRIC whether caching should be symmetric or not
      SSEGeometryPackerCacheWrapper< t_ReturnType>
      (
        const util::ShPtr< SSEGeometryPackerInterface< t_ReturnType> > &SP_PACKER_FUNCTION,
        const bool IS_SYMMETRIC
      ) :
        m_MinimalInterfaceLength( SP_PACKER_FUNCTION->GetMinimalInterfaceLength()),
        m_CacheFunction
        (
          new math::BinaryFunctionCached< SSEGeometryInterface, SSEGeometryInterface, t_ReturnType>
          (
            SP_PACKER_FUNCTION,
            &SSEGeometryInterface::GetGeometryDestructorSignal,
            IS_SYMMETRIC
          )
        )
      {
        // add signal handler for coordinate changes
        m_CacheFunction->AddSignalHandlerForArgument( &SSEGeometryInterface::GetGeometryCoordinateChangeSignal);
      }

      //! @brief Clone function
      //! @return pointer to new SSEGeometryPackerCacheWrapper< t_ReturnType>
      SSEGeometryPackerCacheWrapper< t_ReturnType> *Clone() const
      {
        return new SSEGeometryPackerCacheWrapper< t_ReturnType>( *this);
      }

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
      double GetMinimalInterfaceLength() const
      {
        return m_MinimalInterfaceLength;
      }

    ///////////////
    // operators //
    ///////////////

      //! @brief calculate the packing between all fragments of the given geometries and return it
      //! @param SSE_GEOMETRY_A first SSEGeometryInterface derived class of interest
      //! @param SSE_GEOMETRY_B second SSEGeometryInterface derived class of interest
      //! @return t_ReturnType that contains SSEPackings calculated
      t_ReturnType operator()
      (
        const SSEGeometryInterface &SSE_GEOMETRY_A,
        const SSEGeometryInterface &SSE_GEOMETRY_B
      ) const
      {
        // return from cached function
        return m_CacheFunction->operator()( SSE_GEOMETRY_A, SSE_GEOMETRY_B);
      }

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM)
      {
        // read members
        io::Serialize::Read( m_MinimalInterfaceLength, ISTREAM);
        io::Serialize::Read( m_CacheFunction, ISTREAM);

        // end
        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        // write members
        io::Serialize::Write( m_MinimalInterfaceLength, OSTREAM, INDENT) << '\n';
        io::Serialize::Write( m_CacheFunction, OSTREAM, INDENT);

        // end
        return OSTREAM;
      }

    }; // template class SSEGeometryPackerCacheWrapper

    // instantiate s_Instance
    template< typename t_ReturnType>
    const util::SiPtr< const util::ObjectInterface> SSEGeometryPackerCacheWrapper< t_ReturnType>::s_Instance
    (
      GetObjectInstances().AddInstance( new SSEGeometryPackerCacheWrapper< t_ReturnType>())
    );

  } // namespace assemble
} // namespace bcl

#endif // BCL_ASSEMBLE_SSE_GEOMETRY_PACKER_CACHE_WRAPPER_H_ 
