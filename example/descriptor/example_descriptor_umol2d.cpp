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
#include "descriptor/bcl_descriptor_umol2d.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_fragment_complete.h"
#include "chemistry/bcl_chemistry_molecule_ensemble.h"
#include "io/bcl_io_file.h"
#include "io/bcl_io_ifstream.h"
#include "sdf/bcl_sdf_fragment_factory.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_descriptor_umol2d.cpp
  //!
  //! @author vuot2
  //! @date Feb 8, 2017
  //! @remarks status complete
  //! @remarks reviewed by
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleDescriptorUMol2D :
    public ExampleInterface
  {
  public:

    ExampleDescriptorUMol2D *Clone() const
    {
      return new ExampleDescriptorUMol2D( *this);
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

    // A functor that compare the size of two set of size_t
    class BiggerSet
    {
    public:
      bool operator()
      (
        std::pair< chemistry::AtomEnvironmentBender, storage::Set< size_t> > A,
        std::pair< chemistry::AtomEnvironmentBender, storage::Set< size_t> > B
      )
      {
        return A.second.GetSize() < B.second.GetSize();
      }
    };

    int Run() const
    {
    /////////////////
    // preparation //
    /////////////////

      // Read in small molecule
      io::IFStream input_sdf;
      BCL_ExampleMustOpenInputFile( input_sdf, AddExampleInputPathToFilename( e_Chemistry, "cyclohexane.sdf"));

      // load information into small_mol
      chemistry::FragmentComplete molecule( sdf::FragmentFactory::MakeFragment( input_sdf));

      io::File::CloseClearFStream( input_sdf);

      io::IFStream input;
      BCL_ExampleMustOpenInputFile( input, AddExampleInputPathToFilename( e_Chemistry, "taxol.sdf"));

      // load information into small_mol
      chemistry::FragmentComplete taxol( sdf::FragmentFactory::MakeFragment( input));

      io::File::CloseClearFStream( input);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // test empty constructor
      descriptor::UMol2D molprint2d;

      BCL_ExampleCheck( molprint2d.GetAlias(), "UMol2D");

    /////////////////
    // data access //
    /////////////////

      // check if the default
      const size_t key_size( 574);
      BCL_ExampleCheck( molprint2d.GetNormalSizeOfFeatures(), key_size);

      // Generate the VDW RDF code
      BCL_ExampleCheck( molprint2d( molecule).GetSize(), key_size);
      float expected_vec_array [ key_size] =
      {
          0, 0, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0
      };
      linal::VectorConstReference< float> expected_vector( key_size, expected_vec_array);
      BCL_ExampleCheck( molprint2d( molecule), expected_vector);

      // Generate the VDW RDF code for taxol
      float expected_taxol_array [ key_size] =
      {
         15, 6, 2, 0, 1, 2, 4, 3, 0, 3, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1,
        0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
        0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
        3, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0,
        0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0
      };
      linal::VectorConstReference< float> expected_taxol_vector( key_size, expected_taxol_array);
      BCL_ExampleCheck( molprint2d( taxol), expected_taxol_vector);

    //////////////////////
    // input and output //
    //////////////////////

      // check i/o
      BCL_ExampleIndirectCheck
      (
        TestBCLObjectIOForSymmetry( molprint2d, descriptor::UMol2D()),
        true,
        "UMol2D I/O"
      );

    //////////////////////
    // helper functions //
    //////////////////////

      // This is used for generating list of most common atom environments
      // from a vector of strings of filenames
//      storage::Vector< std::string> filenames;
//
//      filenames.PushBack( AddExampleInputPathToFilename( e_Chemistry, "fixChemCart/ChemCart_mostof.sdf"));
//      filenames.PushBack
//      (
//        AddExampleInputPathToFilename( e_Chemistry, "Vanderbilt_Discovery_Collection_--_Set_1_cleaned.sdf")
//      );
//      filenames.PushBack
//      (
//        AddExampleInputPathToFilename
//        (
//          e_Chemistry,
//          "emolecules_cleaned_3D_1-532_pass_pains_reos_pass_lipinski_unique_labeled_cleaned.sdf"
//        )
//      );
//
//      // creates a fragment feed from filename vector
//      FragmentFeed feed( filenames, sdf::e_Maintain);
//
//      // the map that contains the hashed AEs as keys and the indices of the molecules
//      storage::Map< chemistry::AtomEnvironment2, storage::Set< size_t> > map;
//      size_t molecule_index( 0);
//
//      // read in from the stream until we reach the end of the file or the last index
//      for( ; feed.NotAtEnd(); ++feed)
//      {
//        chemistry::MoleculeEnvironment env_elem( chemistry::AtomEnvironment2::e_AtomRC, *feed);
//        chemistry::MoleculeEnvironment::t_MoleculeEnv mol_env( env_elem.GetMoleculeEnvironment());
//        for( chemistry::MoleculeEnvironment::t_MoleculeEnv::iterator iter = mol_env.Begin(), end_iter = mol_env.End();
//            iter != end_iter; ++iter)
//        {
//          // add the molecule index into the set that corresponds to a certain atom environment key in the map
//          if( map.Has( *iter) == false)
//          {
//            storage::Set< size_t> set;
//            map.Insert( storage::Pair< chemistry::AtomEnvironment2, storage::Set< size_t> >( *iter, set));
//          }
//          map[ *iter].InsertElement( molecule_index);
//        }
//        ++molecule_index;
//      }
//      //std::cout << "There are " << molecule_index << " molecules" << std::endl;
//
//      // store the most common AEs
//      storage::Vector< chemistry::AtomEnvironment2> mostCommonAEs;
//
//      // Store the number of the molecules coverred so far
//      //storage::Vector<size_t> moleculeNum;
//      // sum of numbers of molecules coverred by all the most common AEs so far
//      size_t sum( 0);
//      size_t count( 1);
//
//      // repeat the process as long as the map is not empty
//      while( !map.IsEmpty())
//      {
//
//        // find the most common AEs
//        storage::Map< chemistry::AtomEnvironment2, storage::Set< size_t>>::iterator AE
//            (
//              std::max_element( map.Begin(), map.End(), BiggerSet())
//            );
//        //std::cout << "AE : " << AE->first.UnHash() << std::endl;
//        // add the most common AE to a vector
//        mostCommonAEs.PushBack( AE->first);
//
//        // temporary set contains the molecules covered by AE
//        storage::Set< size_t> tempSet( AE->second);
//
//        // remove the molecule covered by the most common AE from sets of other AEs
//        map.RemoveElement( AE);
//        for( storage::Map< chemistry::AtomEnvironment2, storage::Set< size_t> >::iterator iter = map.Begin();
//            iter != map.End(); ++iter)
//        {
//          storage::Set< size_t> temp;
//          std::set_difference(
//            iter->second.Begin(),
//            iter->second.End(),
//            tempSet.Begin(),
//            tempSet.End(),
//            std::inserter( temp.InternalData(), temp.End()));
//          iter->second.Swap( temp);
//        }
//        // store the number of molecules which are covered by the most common AEs
//        sum += tempSet.GetSize();
//        //std::cout << count << " sum = " << sum << std::endl;
//        std::cout << count << " sum = " << ( sum * 1.0 / molecule_index * 100.0) << std::endl;
//        ++count;
//        moleculeNum.PushBack(sum/molecule_index * 100);
//      }

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleDescriptorMolprint2D

  const ExampleClass::EnumType ExampleDescriptorUMol2D::s_Instance
  (
    GetExamples().AddEnum( ExampleDescriptorUMol2D())
  );

} // namespace bcl

