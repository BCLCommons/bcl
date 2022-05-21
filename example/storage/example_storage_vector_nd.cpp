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
#include "storage/bcl_storage_vector_nd.h"

// includes from bcl - sorted alphabetically
#include "coord/bcl_coord_axes.h"
#include "util/bcl_util_si_ptr_vector.h"
#include "util/bcl_util_stopwatch.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_storage_vector_nd.cpp
  //!
  //! @author alexanns
  //! @date
  //! @remarks status empty
  //! @remarks reviewed by nobody on 
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleStorageVectorND :
    public ExampleInterface
  {
  public:

    ExampleStorageVectorND *Clone() const
    { return new ExampleStorageVectorND( *this);}

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
      //constructor and access function for const containers
      {
        linal::Vector3D a( 1.1), b( 2.2), c( 3.3), d( 4.4);
        const storage::VectorND< 1, util::SiPtr< const linal::Vector3D> > container1( a);
        const storage::VectorND< 2, util::SiPtr< const linal::Vector3D> > container2( a, b);
        const storage::VectorND< 3, util::SiPtr< const linal::Vector3D> > container3( a, b, c);
        const storage::VectorND< 4, util::SiPtr< const linal::Vector3D> > container4( a, b, c, d);

        BCL_MessageStd( "use access function First(), Second(), Third() and Fourth()");
        BCL_MessageStd( "const container with one element:    |" + util::Format()( container1.First()) + "|");
        BCL_MessageStd( "const container with two elements:   |" + util::Format()( container2.First()) + "|" + util::Format()( container2.Second()) + "|");
        BCL_MessageStd( "const container with three elements: |" + util::Format()( container3.First()) + "|" + util::Format()( container3.Second()) + "|" + util::Format()( container3.Third()) + "|");
        BCL_MessageStd( "const container with four elements:  |" + util::Format()( container4.First()) + "|" + util::Format()( container4.Second()) + "|" + util::Format()( container4.Third()) + "|" + util::Format()( container4.Fourth()) + "|");

        BCL_MessageStd( "use access function operator()( POS)");
        BCL_MessageStd( "const container with one element:    |" + util::Format()( container1( 0)) + "|");
        BCL_MessageStd( "const container with two elements:   |" + util::Format()( container2( 0)) + "|" + util::Format()( container2( 1)) + "|");
        BCL_MessageStd( "const container with three elements: |" + util::Format()( container3( 0)) + "|" + util::Format()( container3( 1)) + "|" + util::Format()( container3( 2)) + "|");
        BCL_MessageStd( "const container with four elements:  |" + util::Format()( container4( 0)) + "|" + util::Format()( container4( 1)) + "|" + util::Format()( container4( 2)) + "|" + util::Format()( container4( 3)) + "|");
      }

      //constructor and access function for containers
      {
        storage::VectorND< 1, double> container1( 1.1);
        storage::VectorND< 2, double> container2( 1.1, 2.2);
        storage::VectorND< 3, double> container3( 1.1, 2.2, 3.3);
        storage::VectorND< 4, double> container4( 1.1, 2.2, 3.3, 4.4);

        BCL_MessageStd( "use access function First(), Second(), Third() and Fourth()");
        BCL_MessageStd( "container with one element:    |" + util::Format()( container1.First()) + "|");
        BCL_MessageStd( "container with two elements:   |" + util::Format()( container2.First()) + "|" + util::Format()( container2.Second()) + "|");
        BCL_MessageStd( "container with three elements: |" + util::Format()( container3.First()) + "|" + util::Format()( container3.Second()) + "|" + util::Format()( container3.Third()) + "|");
        BCL_MessageStd( "container with four elements:  |" + util::Format()( container4.First()) + "|" + util::Format()( container4.Second()) + "|" + util::Format()( container4.Third()) + "|" + util::Format()( container4.Fourth()) + "|");

        BCL_MessageStd( "use access function operator()( POS)");
        BCL_MessageStd( "container with one element:    |" + util::Format()( container1( 0)) + "|");
        BCL_MessageStd( "container with two elements:   |" + util::Format()( container2( 0)) + "|" + util::Format()( container2( 1)) + "|");
        BCL_MessageStd( "container with three elements: |" + util::Format()( container3( 0)) + "|" + util::Format()( container3( 1)) + "|" + util::Format()( container3( 2)) + "|");
        BCL_MessageStd( "container with four elements:  |" + util::Format()( container4( 0)) + "|" + util::Format()( container4( 1)) + "|" + util::Format()( container4( 2)) + "|" + util::Format()( container4( 3)) + "|");

        BCL_MessageStd( "change last element in each container: ");
        container1( 0) = 11.11;
        BCL_MessageStd( util::Format()( container1( 0)));
        container2( 1) = 22.22;
        BCL_MessageStd( util::Format()( container2( 1)));
        container3( 2) = 33.33;
        BCL_MessageStd( util::Format()( container3( 2)));
        container4( 3) = 44.44;
        BCL_MessageStd( util::Format()( container4( 3)));
      }

      //compare StorageVectorND and util::SiPtrVector in time performance
      //construct StorageVectorND
      BCL_MessageStd( "compare util::SiPtrVector to storage::VectorND< 2, util::SiPtr> for 10000 constructions: ");
      {
        util::Stopwatch stopwatch;
        util::SiPtr< const util::Stopwatch> sp( &stopwatch);
        for( size_t i( 0); i < 100000; ++i)
        {
          util::SiPtrVector< const util::Stopwatch>( sp, sp);
        }
      }

      {
        util::Stopwatch stopwatch;
        util::SiPtr< const util::Stopwatch> sp( &stopwatch);
        for( size_t i( 0); i < 100000; ++i)
        {
          storage::VectorND< 2, util::SiPtr< const util::Stopwatch> >( sp, sp);
        }
      }

      return 0;
    } //end ExampleClass::ExampleStorageVectorND

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleStorageVectorND

  const ExampleClass::EnumType ExampleStorageVectorND::s_Instance
  (
    GetExamples().AddEnum( ExampleStorageVectorND())
  );

} // namespace bcl
