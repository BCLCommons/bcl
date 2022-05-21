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

#ifndef BCL_UTIL_BINARY_FUNCTION_STL_WRAPPER_H_
#define BCL_UTIL_BINARY_FUNCTION_STL_WRAPPER_H_

// include the namespace header
#include "bcl_util.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_util_binary_function_interface_serializable.h"
#include "bcl_util_object_instances.h"
#include "io/bcl_io_serialization.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace util
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class BinaryFunctionSTLWrapper
    //! @brief wraps stl functions into a bcl binary function interface
    //! @tparam t_BinaryFunction a binary function such as std::less, std::equal_to, etc.
    //!
    //! @see @link example_util_binary_function_stl_wrapper.cpp @endlink
    //! @author woetzen
    //! @date 01.04.2007
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_BinaryFunction>
    class BinaryFunctionSTLWrapper :
      public BinaryFunctionInterfaceSerializable
             <
               typename t_BinaryFunction::first_argument_type,
               typename t_BinaryFunction::second_argument_type,
               typename t_BinaryFunction::result_type
             >
    {

    private:

    //////////
    // data //
    //////////

      t_BinaryFunction m_Function;

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const SiPtr< const ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! default constructor
      BinaryFunctionSTLWrapper() :
        m_Function()
      {
      }

      //! @brief construct from a binary function
      //! @param BINARY_FUNCTION binary function to be used instead of default constructed
      BinaryFunctionSTLWrapper( const t_BinaryFunction &BINARY_FUNCTION) :
        m_Function( BINARY_FUNCTION)
      {
      }

      //! virtual copy constructor
      BinaryFunctionSTLWrapper *Clone() const
      {
        return new BinaryFunctionSTLWrapper( *this);
      }

      //! @brief operator the implements the binary operation on the two arguments returning a result
      //! @param ARGUMENT1 argument 1
      //! @param ARGUMENT2 argument 2
      //! @return the Result of the binary operation
      typename t_BinaryFunction::result_type operator()
      (
        const typename t_BinaryFunction::first_argument_type &ARGUMENT1,
        const typename t_BinaryFunction::second_argument_type &ARGUMENT2
      ) const
      {
        return m_Function( ARGUMENT1, ARGUMENT2);
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

      //! @brief returns parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      virtual io::Serializer GetSerializer() const
      {
        io::Serializer serializer;
        serializer.SetClassDescription( "This class needs a proper implementation of this function.");
        return serializer;
      }

    //////////////////////
    // input and output //
    //////////////////////

    }; // template class BinaryFunctionSTLWrapper

    // instantiate s_Instance
    template< typename t_BinaryFunction>
    const SiPtr< const ObjectInterface> BinaryFunctionSTLWrapper< t_BinaryFunction>::s_Instance
    (
      GetObjectInstances().AddInstance( new BinaryFunctionSTLWrapper< t_BinaryFunction>())
    );

  } // namespace util
} // namespace bcl

#endif // BCL_UTIL_BINARY_FUNCTION_STL_WRAPPER_H_ 
