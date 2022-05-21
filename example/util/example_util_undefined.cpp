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
#include "util/bcl_util_undefined.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_util_undefined.cpp
  //!
  //! @author woetzen
  //! @date   Nov 06, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleUtilUndefined :
    public ExampleInterface
  {
  public:

    ExampleUtilUndefined *Clone() const
    { return new ExampleUtilUndefined( *this);}

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    int Run() const
    {
    ///////////////////////
    // undefined objects //
    ///////////////////////

      BCL_MessageStd( "double Undefined: " + util::Format()( util::GetUndefined< double>()));
      BCL_MessageStd( "double Undefined: " + util::Format()( util::GetUndefined< double>()));
      BCL_MessageStd( "float  Undefined: " + util::Format()( util::GetUndefined<  float>()));
      BCL_MessageStd( "int    Undefined: " + util::Format()( util::GetUndefined<    int>()));
      BCL_MessageStd( "size_t Undefined: " + util::Format()( util::GetUndefined< size_t>()));
      BCL_MessageStd( "long   Undefined: " + util::Format()( util::GetUndefined<   long>()));
      BCL_MessageStd( "short  Undefined: " + util::Format()( util::GetUndefined<  short>()));

      //double
      BCL_MessageStd
      (
        "is util::GetUndefined< double>() defined?: " + util::Format()( util::IsDefined( util::GetUndefined< double>()))
      );
      BCL_Example_Check
      (
        !util::IsDefined( util::GetUndefined< double>()),
        "util::GetUndefined< double>() has to return false in function call IsDefined()"
      );

      //float
      BCL_MessageStd
      (
        "is util::GetUndefined< float>() defined?: " + util::Format()( util::IsDefined( util::GetUndefined< float>()))
      );
      BCL_Example_Check
      (
        !util::IsDefined( util::GetUndefined< float>()),
        "util::GetUndefined< float>() has to return false in function call IsDefined()"
      );

      //int
      BCL_MessageStd
      (
        "is util::GetUndefined<    int>() defined?: " + util::Format()( util::IsDefined( util::GetUndefined<    int>()))
      );
      BCL_Example_Check
      (
        !util::IsDefined( util::GetUndefined<   int>()),
        "util::GetUndefined<    int>() has to return false in function call IsDefined()"
      );

      //size_t
      BCL_MessageStd
      (
        "is util::GetUndefined< size_t>() defined?: " + util::Format()( util::IsDefined( util::GetUndefined< size_t>()))
      );
      BCL_Example_Check
      (
        !util::IsDefined( util::GetUndefined< size_t>()),
        "util::GetUndefined< size_t>() has to return false in function call IsDefined()"
      );

      //long
      BCL_MessageStd
      (
        "is util::GetUndefined<   long>() defined?: " + util::Format()( util::IsDefined( util::GetUndefined<   long>()))
      );
      BCL_Example_Check
      (
        !util::IsDefined( util::GetUndefined<   long>()),
        "util::GetUndefined<   long>() has to return false in function call IsDefined()"
      );

      //short
      BCL_MessageStd
      (
        "is util::GetUndefined<  short>() defined?: " + util::Format()( util::IsDefined( util::GetUndefined< short>()))
      );
      BCL_Example_Check
      (
        !util::IsDefined( util::GetUndefined< short>()),
        "util::GetUndefined<  short>() has to return false in function call IsDefined()"
      );

    /////////////////////
    // defined objects //
    /////////////////////

      //double
      BCL_MessageStd
      (
        "is double 5.5 defined?: " + util::Format()( util::IsDefined( double( 5.5)))
      );
      BCL_Example_Check
      (
        util::IsDefined( double( 5.5)), "double 5.5 is not defined. Has to return true for call IsDefined( double)"
      );

      //float
      BCL_MessageStd
      (
        "is  float 6.6 defined?: " + util::Format()( util::IsDefined( float( 6.6)))
      );
      BCL_Example_Check
      (
        util::IsDefined( float( 6.6)), " float 6.6 is not defined. Has to return true for call IsDefined(  float)"
      );

      //int
      BCL_MessageStd
      (
        "is    int  -6 defined?: " + util::Format()( util::IsDefined(    int(  -6)))
      );
      BCL_Example_Check
      (
        util::IsDefined(    int(  -6)), "   int  -6 is not defined. Has to return true for call IsDefined(    int)"
      );

      //size_t
      BCL_MessageStd
      (
        "is size_t   7 defined?: " + util::Format()( util::IsDefined( size_t(   7)))
      );
      BCL_Example_Check
      (
        util::IsDefined( size_t(   7)), "size_t   7 is not defined. Has to return true for call IsDefined( size_t)"
      );

      //long
      BCL_MessageStd
      (
        "is   long  77 defined?: " + util::Format()( util::IsDefined(   long(  77)))
      );
      BCL_Example_Check
      (
        util::IsDefined(   long(  77)), "long    77 is not defined. Has to return true for call IsDefined(   long)"
      );

      //short
      BCL_MessageStd
      (
        "is   short  5 defined?: " + util::Format()( util::IsDefined(  short(   5))) + "\n"
      );
      BCL_Example_Check
      (
        util::IsDefined(  short(   7)), "short    7 is not defined. Has to return true for call IsDefined(   short)"
      );

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleUtilUndefined

  const ExampleClass::EnumType ExampleUtilUndefined::s_Instance
  (
    GetExamples().AddEnum( ExampleUtilUndefined())
  );

} // namespace bcl
