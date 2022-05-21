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

#ifndef BCL_UTIL_WRAPPER_BASE_H_
#define BCL_UTIL_WRAPPER_BASE_H_

// include the namespace header
#include "bcl_util.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_util_object_instances.h"
#include "bcl_util_object_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace util
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class WrapperBase
    //! @brief This is a template class wrapping base types like int, double.. without class behavior, so that they have a
    //! common bcl interface.
    //! @details They will still have all functionalities, since there is an operator converting it to its DATA_TYPE.
    //!
    //! @see @link example_util_wrapper_base.cpp @endlink
    //! @author woetzen
    //! @date 29.05.2007
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_BaseType>
    class WrapperBase :
      public ObjectInterface
    {

    private:

    //////////
    // data //
    //////////

      t_BaseType m_Data;

    public:

      //! single instance of that class
      static const SiPtr< const ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      WrapperBase< t_BaseType>() :
        m_Data()
      {
      }

      //! @brief constructor from DATA
      //! @param DATA ref to instance of t_BaseType
      WrapperBase< t_BaseType>( const t_BaseType &DATA) :
        m_Data( DATA)
      {
      }

      //! @brief copy constructor
      //! @param WRAPPER wrapper object to be copied
      WrapperBase< t_BaseType>( const WrapperBase< t_BaseType> &WRAPPER) :
        m_Data( WRAPPER)
      {
      }

      //! @brief virtual copy constructor
      WrapperBase< t_BaseType> *Clone() const
      {
        return new WrapperBase< t_BaseType>( *this);
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

    ///////////////
    // operators //
    ///////////////

      //! @brief for implicit conversion into the base type
      operator t_BaseType()
      {
        return m_Data;
      }

      //! @brief for implicit conversion into const base type
      operator const t_BaseType() const
      {
        return m_Data;
      }

    //////////////////////
    // input and output //
    //////////////////////
    protected:

      //! read from std::istream
      std::istream &Read( std::istream &ISTREAM)
      {
        return ISTREAM >> m_Data;
      }

      //! write to std::ostream
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        return OSTREAM << m_Data;
      }

    }; // class WrapperBase

    template< typename t_ArgumentType>
    const SiPtr< const ObjectInterface> WrapperBase< t_ArgumentType>::s_Instance
    (
      GetObjectInstances().AddInstance( new WrapperBase< t_ArgumentType>())
    );

  } // namespace util
} // namespace bcl

#endif //BCL_UTIL_WRAPPER_BASE_H_
