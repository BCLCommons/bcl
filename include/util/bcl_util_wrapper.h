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

#ifndef BCL_UTIL_WRAPPER_H_
#define BCL_UTIL_WRAPPER_H_

// include the namespace header
#include "bcl_util.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_util_class_descriptor.h"
#include "bcl_util_object_instances.h"
#include "bcl_util_object_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace util
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Wrapper
    //! @brief class wrapping non bcl classes into for use with the common bcl interface
    //! @details This is a template class wrapping non bcl classes, so that they have a common bcl interface.
    //! They will still have all functionalities, since the class is derived from the wrapped class.
    //!
    //! @see @link example_util_wrapper.cpp @endlink
    //! @author woetzen
    //! @date 29.05.2007
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DataType>
    class Wrapper :
      public ObjectInterface,
      public t_DataType
    {

    private:

    //////////
    // data //
    //////////

    public:

      //! single instance of that class
      static const SiPtr< const ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //!default constructor
      Wrapper() :
        t_DataType()
      {
      }

      //!constructor from t_DataType
      Wrapper( const t_DataType &DATA) :
        t_DataType( DATA)
      {
      }

      //!copy constructor
      Wrapper( const Wrapper< t_DataType> &WRAPPER) :
        t_DataType( WRAPPER)
      {
      }

      //!Clone
      Wrapper< t_DataType> *Clone() const
      {
        return new Wrapper< t_DataType>( *this);
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

      //! convert this to a reference of t_DataType
      t_DataType &GetData()
      {
        return *this;
      }

      //! convert this to a const reference of t_DataType
      t_DataType const &GetData() const
      {
        return *this;
      }

    //////////////////////
    // input and output //
    //////////////////////
    protected:

      //! read from std::istream
      std::istream &Read( std::istream &ISTREAM)
      {
        return ISTREAM >> GetData();
      }

      //! write to std::ostream
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        return OSTREAM << GetData();
      }

    }; // class Wrapper

    template< typename t_ArgumentType>
    const SiPtr< const ObjectInterface> Wrapper< t_ArgumentType>::s_Instance
    (
      GetObjectInstances().AddInstance( new Wrapper< t_ArgumentType>())
    );

  } // namespace util
} // namespace bcl

#endif //BCL_UTIL_WRAPPER_H_
