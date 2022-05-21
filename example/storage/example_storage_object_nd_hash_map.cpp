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
#include "storage/bcl_storage_object_nd_hash_map.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_storage_object_nd_hash_map.cpp
  //!
  //! @author karakam
  //! @date
  //! @remarks status empty
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleStorageObjectNDHashMap :
    public ExampleInterface
  {
  public:

    ExampleStorageObjectNDHashMap *Clone() const
    { return new ExampleStorageObjectNDHashMap( *this);}

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
      // test default constructor
      storage::ObjectNDHashMap< 2, storage::Vector< double>, double> map;

      // test Insert function
      storage::Vector< double> sv1( 4, 0.2);
      storage::Vector< double> sv2( 5, 0.1);
      map.Insert( storage::VectorND< 2, util::SiPtr< const storage::Vector< double> > >( sv1, sv2), 0.4);
      map.Insert( storage::VectorND< 2, util::SiPtr< const storage::Vector< double> > >( sv1, sv2), 0.5);
      for( storage::HashMap< size_t, double>::const_iterator itr( map.Begin()), itr_end( map.End()); itr != itr_end; ++itr)
      {
        BCL_MessageStd( util::Format()( itr->second));
      }

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleStorageObjectNDHashMap

  const ExampleClass::EnumType ExampleStorageObjectNDHashMap::s_Instance
  (
    GetExamples().AddEnum( ExampleStorageObjectNDHashMap())
  );

} // namespace bcl
