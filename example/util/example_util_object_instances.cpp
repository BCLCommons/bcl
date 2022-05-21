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
#include "util/bcl_util_object_instances.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_file.h"
#include "math/bcl_math_template_instantiations.h"
#include "math/bcl_math_tensor.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //
  //! @example example_util_object_instances.cpp
  //!
  //! @author woetzen
  //! @date Nov 06, 2010
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleUtilObjectInstances :
    public ExampleInterface
  {
  public:

    ExampleUtilObjectInstances *Clone() const
    { return new ExampleUtilObjectInstances( *this);}

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
//      for
//      (
//        util::ObjectInstances::const_iterator
//          itr( GetObjectInstances().Begin()), itr_end( GetObjectInstances().End());
//        itr!= itr_end; ++itr)
//      {
//        util::GetLogger()<< itr->GetName() << '\n';
//      }

      BCL_MessageStd
      (
        "number of object instances available: " +
        util::Format()( GetObjectInstances().GetKnownObjectNames().GetSize())
      );

      //write serialized objects
      util::ShPtr< linal::VectorInterface< double> >       vector( new linal::Vector< double>( 6, 4.0));
      util::ShPtr< linal::MatrixInterface< double> >       matrix( new linal::Matrix< double>( 2, 3, 4.0));
      util::ShPtr< math::Tensor< double> >                 tensor( new math::Tensor< double>( 2, 2, 3, 4.0));

      io::OFStream write;
      std::string filename( AddExampleOutputPathToFilename( util::GetNamespaceIdentifier(), "serialization.test"));
      BCL_ExampleMustOpenOutputFile( write, filename);

      io::Serialize::Write( vector, write) << '\n';
      io::Serialize::Write( matrix, write) << '\n';
      io::Serialize::Write( tensor, write);

      io::File::CloseClearFStream( write);

      //read serialized objects
      util::Format format;
      io::IFStream read;
      BCL_ExampleMustOpenInputFile( read, filename);

      util::ShPtr< linal::VectorInterface< double> > object1;
      read >> object1;
      BCL_MessageStd( format( object1));
      util::ShPtr< linal::MatrixInterface< double> > object2;
      read >> object2;
      BCL_MessageStd( format( object2));
      util::ShPtr< math::Tensor< double> > object3;
      read >> object3;
      BCL_MessageStd( format( object3));

      io::File::CloseClearFStream( read);

      return 0;
    } //end ExampleClass::ExampleUtilObjectInstances

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleUtilObjectInstances

  const ExampleClass::EnumType ExampleUtilObjectInstances::s_Instance
  (
    GetExamples().AddEnum( ExampleUtilObjectInstances())
  );

} // namespace bcl
