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

#ifndef BCL_MATH_IDENTITY_H_
#define BCL_MATH_IDENTITY_H_

// include the namespace header
#include "bcl_math.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_math_function_interface_serializable.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Identity
    //! @brief This class is derived from function interface and is designed to give back
    //!        the same value as passed to the operator also called 'identity'.
    //!
    //! @tparam t_DataType template parameter for data type of passed argument
    //!
    //! @see @link example_math_identity.cpp @endlink
    //! @author butkiem1
    //! @date Dec 2, 2009
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DataType>
    class Identity :
      public FunctionInterfaceSerializable< t_DataType, t_DataType>
    {

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief Clone function
      //! @return pointer to new Identity< t_DataType>
      Identity< t_DataType> *Clone() const
      {
        return new Identity< t_DataType>( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      //! @brief gets the scheme for this function class
      //! @return the scheme for this function class
      const std::string &GetScheme() const
      {
        static const std::string s_scheme( "identity");
        return s_scheme;
      }

    ///////////////
    // operators //
    ///////////////

      //! @brief operator taking an ARGUMENT and returning it as the result
      //! @param ARGUMENT Argument to be used to evaluate the function
      //! @return function value of the given argument
      t_DataType operator()( const t_DataType &ARGUMENT) const
      {
        return ARGUMENT;
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
        // end
        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        // end
        return OSTREAM;
      }

    }; // template class Identity

    // instantiate s_Instance
    template< typename t_DataType>
    const util::SiPtr< const util::ObjectInterface> Identity< t_DataType>::s_Instance
    (
      GetObjectInstances().AddInstance( new Identity< t_DataType>())
    );

  } // namespace math
} // namespace bcl

#endif // BCL_MATH_IDENTITY_H_ 
