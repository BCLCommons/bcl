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

// include example header
#include "example.h"
// include the header of the class which this example is for
#include "type/bcl_type_enable_if.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_type_enable_if.cpp
  //!
  //! @author mendenjl
  //! @date Nov 29, 2012
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleTypeEnableIf :
    public ExampleInterface
  {
  public:

    ExampleTypeEnableIf *Clone() const
    { return new ExampleTypeEnableIf( *this);}

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief example of a function using EnableIf as a return type
    //! @return true; this function will be called if t_DataType is built-in
    template< typename t_DataType>
    typename type::EnableIf< std::numeric_limits< t_DataType>::is_specialized, bool>::Type
      IsIntrinsicType( t_DataType) const
    {
      return true;
    }

    //! @brief example of a function using EnableIf as a return type
    //! @return false; this function will be called if t_DataType is not built-in
    template< typename t_DataType>
    typename type::EnableIf< !std::numeric_limits< t_DataType>::is_specialized, bool>::Type
      IsIntrinsicType( t_DataType) const
    {
      return false;
    }

    int Run() const
    {
    //////////////
    // EnableIf //
    //////////////

      BCL_ExampleCheck( IsIntrinsicType( int( 1)), true);
      BCL_ExampleCheck( IsIntrinsicType( std::vector< float>()), false);

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleTypeEnableIf

  const ExampleClass::EnumType ExampleTypeEnableIf::s_Instance
  (
    GetExamples().AddEnum( ExampleTypeEnableIf())
  );

} // namespace bcl
