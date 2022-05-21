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
#include "util/bcl_util_class_descriptor.h"
#include "util/bcl_util_class_name_standardizer.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

class classbase
{
};

namespace bcl
{
  class Base
  {
  public:

    //! destructor
    virtual ~Base()
    {
    }

    //! @brief returns class name
    //! @return the class name as const ref std::string
    virtual const std::string &GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }
  };

  class Derived :
    public Base
  {
  public:

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }
  };

  class DerivedV :
    public virtual Base
  {
  public:

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }
  };

  class DerivedV2 :
    public virtual Base,
    public virtual DerivedV
  {
  public:

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }
  };

  template< typename t_Data, typename t_MoreData>
  class DerivedV3 :
    public virtual DerivedV,
    public virtual DerivedV2
  {
  public:

    enum EnumTest
    {
      yes = 1,
      no  = 0
    };

    union UnionTest
    {
      double a;
      int b;
    };

    struct StructTest
    {
    };

    t_Data     *m_Unused;
    t_MoreData m_Unused2;

  public:

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }
  };

  template< DerivedV3< double, double>::EnumTest YN>
  class EnumTemplateClass
  {
    enum
    {
      MyValue = YN
    };

  };

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //
  //! @example example_util_class_descriptor.cpp
  //! @details tests that GetStaticClassName returns the expected names for classes of large complexity
  //!
  //! @author mendenjl
  //! @date July 01, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleUtilClassDescriptor :
    public ExampleInterface
  {
  public:

    ExampleUtilClassDescriptor *Clone() const
    {
      return new ExampleUtilClassDescriptor( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &GetClassIdentifier() const
    {
      return GetStaticClassName< ExampleUtilClassDescriptor>();
    }

    int Run() const
    {
    //////////////////////
    // helper functions //
    //////////////////////

      // test GetStaticClassName on all basic types
      BCL_ExampleCheck( GetStaticClassName< long double>(), "long-double");
      BCL_ExampleCheck( GetStaticClassName< double>(), "double");
      BCL_ExampleCheck( GetStaticClassName< float>(), "float");
      BCL_ExampleCheck( GetStaticClassName< void>(), "void");
      BCL_ExampleCheck( GetStaticClassName< bool>(), "bool");
      BCL_ExampleCheck( GetStaticClassName< short>(), "short");
      BCL_ExampleCheck( GetStaticClassName< int>(), "int");
      BCL_ExampleCheck( GetStaticClassName< long>(), "long");
      BCL_ExampleCheck( GetStaticClassName< long long>(), "long-long");
      BCL_ExampleCheck( GetStaticClassName< wchar_t>(), "wchar_t");
      BCL_ExampleCheck( GetStaticClassName< char>(), "char");
      BCL_ExampleCheck( GetStaticClassName< size_t>(), "size_t");
      BCL_ExampleCheck( GetStaticClassName< std::string>(), "std::string");

      // see what happens when we make add pointers
      BCL_ExampleCheck( GetStaticClassName< double *>(), "double*");
      BCL_ExampleCheck( GetStaticClassName< int *>(), "int*");
      BCL_ExampleCheck( GetStaticClassName< void *>(), "void*");

      // see what happens when we make add references
      BCL_ExampleCheck( GetStaticClassName< double &>(), "double&");
      BCL_ExampleCheck( GetStaticClassName< int *&>(), "int*&");
      BCL_ExampleCheck( GetStaticClassName< void * &>(), "void*&");

      // test GetStaticClassName usage with const
      BCL_ExampleCheck
      (
        GetStaticClassName< const double *>(),
        "const-double*"
      );
      BCL_ExampleCheck
      (
        GetStaticClassName< const double *const>(),
        "const-double*const"
      );
      BCL_ExampleCheck
      (
        GetStaticClassName< const size_t>(),
        "const-size_t"
      );
      BCL_ExampleCheck
      (
        GetStaticClassName< const void *>(),
        "const-void*"
      );
      BCL_ExampleCheck
      (
        GetStaticClassName< const std::string>(),
        "const-std::string"
      );
      BCL_ExampleCheck
      (
        GetStaticClassName< const unsigned char>(),
        "const-unsigned-char"
      );
      BCL_ExampleCheck
      (
        GetStaticClassName< const unsigned long long>(),
        "const-" + GetStaticClassName< unsigned long long>()
      );
      BCL_ExampleCheck
      (
        GetStaticClassName< const unsigned long long *>(),
        "const-" + GetStaticClassName< unsigned long long>() + "*"
      );

      // test that volatile works, even with const and pointers
      BCL_ExampleCheck
      (
        GetStaticClassName< volatile long>(),
        "volatile-long"
      );
      BCL_ExampleCheck
      (
        GetStaticClassName< const volatile long>(), "const-volatile-long"
      );
      BCL_ExampleCheck
      (
        GetStaticClassName< const volatile long *>(), "const-volatile-long*"
      );
      BCL_ExampleCheck
      (
        GetStaticClassName< volatile long *>(), "volatile-long*"
      );
      BCL_ExampleCheck
      (
        GetStaticClassName< const double *const *volatile const>(),
        "const-double*const*const-volatile"
      );
      BCL_ExampleCheck
      (
        GetStaticClassName< const double *const *const volatile>(),
        "const-double*const*const-volatile"
      );
      BCL_ExampleCheck
      (
        GetStaticClassName< const volatile unsigned long long *const *const volatile>(),
        "const-volatile-" + GetStaticClassName< unsigned long long>() + "*const*const-volatile"
      );

      // Check that get static class name works with regular classes as well
      BCL_ExampleCheck( GetStaticClassName< const Base>(), "const-bcl::Base");
      BCL_ExampleCheck( GetStaticClassName< Base *>(), "bcl::Base*");
      BCL_ExampleCheck
      (
        GetStaticClassName< const volatile Base *const>(),
        "const-volatile-bcl::Base*const"
      );

      // In MSVS 2008, __PRETTY_LINE__ would return something like "classclassbase*", because MSVC omits the space after
      // class,struct,enum, or union if it is a pointer to type (this is probably a bug)
      // Check that our code handles this case correctly
      BCL_ExampleCheck( GetStaticClassName< classbase *>(), "classbase*");

      // Check that templates are handled properly, and that std::string is not expanded to std::basic_string...
      // Due to macro expansion rules, we need to surround arguments-that-have-commas in macros with ()
      BCL_ExampleCheck
      (
        ( GetStaticClassName< DerivedV3< std::string, DerivedV3< const double, std::string> > >()),
        "bcl::DerivedV3<std::string,bcl::DerivedV3<const-double,std::string>>"
      );

      // Ensure that nested classes, enums, and structs are handled properly
      BCL_ExampleCheck
      (
        ( GetStaticClassName< DerivedV3< std::string, DerivedV3< const double, std::string> >::EnumTest>()),
        "bcl::DerivedV3<std::string,bcl::DerivedV3<const-double,std::string>>::EnumTest"
      );
      BCL_ExampleCheck
      (
        ( GetStaticClassName< DerivedV3< std::string, DerivedV3< const double, std::string> >::UnionTest>()),
        "bcl::DerivedV3<std::string,bcl::DerivedV3<const-double,std::string>>::UnionTest"
      );
      BCL_ExampleCheck
      (
        ( GetStaticClassName< DerivedV3< std::string, DerivedV3< const double, std::string> >::StructTest>()),
        "bcl::DerivedV3<std::string,bcl::DerivedV3<const-double,std::string>>::StructTest"
      );
      BCL_ExampleCheck
      (
        ( GetStaticClassName< DerivedV3< std::string, DerivedV3< const double, std::string> >::StructTest *>()),
        "bcl::DerivedV3<std::string,bcl::DerivedV3<const-double,std::string>>::StructTest*"
      );

      // test util::StandardizeClassName
      // check that const is always put before the type (or the nearest pointer or reference, whichever comes first)
      BCL_ExampleCheck
      (
        util::StandardizeClassName( "const double *"),
        util::StandardizeClassName( "double const*")
      );

      // make sure that const is reordered properly even if it is formatted in windows style
      BCL_ExampleCheck
      (
        util::StandardizeClassName( "test<unsigned short const *,otest< double ** const> > const * &const"),
        "const-test<const-unsigned-short*,otest<double**const>>*const&"
      );

      // Templating on enum types is not currently done in the bcl.  When it is, it will become important to test this
      //BCL_ExampleCheck
      //(
      //  (GetStaticClassName< EnumTemplateClass< DerivedV3< double, double>::yes> >()),
      //  "bcl::EnumTemplateClass<(bcl::DerivedV3::EnumTest)1>",
      //  ""
      //);

      // test whether enum types are abbreviated properly
      BCL_ExampleCheck
      (
        util::StandardizeClassName
        (
          "util::Enum<SomeLong<Possibly,Templated>,util::Typethat<may<have<nested<templates>>>>>-const"
        ),
        "const-util::Typethat<may<have<nested<templates>>>>::Enum"
      );

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleUtilClassDescriptor

  const ExampleClass::EnumType
    ExampleUtilClassDescriptor::s_Instance( GetExamples().AddEnum( ExampleUtilClassDescriptor()));
} // namespace bcl

