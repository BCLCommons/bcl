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

//// (c) Copyright BCL @ Vanderbilt University 2014
//// (c) BCL Homepage: http://www.meilerlab.org/bclcommons
//// (c) The BCL software is developed by the contributing members of the BCL @ Vanderbilt University
//// (c) This file is part of the BCL software suite and is made available under license.
//// (c) To view or modify this file, you must enter into one of the following agreements if you have not done so already:
//// (c) For academic and non-profit users:
//// (c)   the BCL Academic Single-User License, available at http://www.meilerlab.org/bclcommons/license
//// (c) For commercial users:
//// (c)   The BCL Commercial Site License, available upon request from bcl-support-commercial@meilerlab.org
//// (c) For BCL developers at Vanderbilt University:
//// (c)   The BCL Developer Agreement, available at http://www.meilerlab.org/bclcommons/developer_agreement
//// (c)
//// (c)   As part of all such agreements, this copyright notice must appear, verbatim and without addition, at the
//// (c) top of all source files of the BCL project and may not be modified by any party except the BCL developers at
//// (c) Vanderbilt University.
//// (c)   The BCL copyright and license yields to non-BCL copyrights and licenses where indicated by code comments.
//// (c)   Questions about this copyright notice or license agreement may be emailed to bcl-support-academic@meilerlab.org
//// (c) (for academic users) or bcl-support-commercial@meilerlab.org (for commercial users)
//
//// include example header
//#include "example.h"
//
//// include the header of the class which this example is for
//#include "chemistry/bcl_chemistry_fragment_remove_bond.h"
//
//// includes from bcl - sorted alphabetically
//#include "chemistry/bcl_chemistry_atom_conformational_interface.h"
//#include "chemistry/bcl_chemistry_collector_valence.h"
//#include "chemistry/bcl_chemistry_pick_atom_random.h"
//#include "chemistry/bcl_chemistry_pick_fragment_random.h"
//#include "io/bcl_io_file.h"
//#include "random/bcl_random_uniform_distribution.h"
//
//// external includes - sorted alphabetically
//#include <stdio.h>
//namespace bcl
//{
//  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  //!
//  //! @example example_chemistry_fragment_remove_bond.cpp
//  //!
//  //! @author brownbp1, mendenjl
//  //! @date Sep 14, 2019
//  //! @remarks status complete
//  //! @remarks reviewed by nobody on
//  //!
//  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  class ExampleChemistryFragmentRemoveBond :
//    public ExampleInterface
//  {
//  public:
//
//    ExampleChemistryFragmentRemoveBond *Clone() const
//    {
//      return new ExampleChemistryFragmentRemoveBond( *this);
//    }
//
//  /////////////////
//  // data access //
//  /////////////////
//
//    //! @brief returns class name
//    //! @return the class name as const ref std::string
//    const std::string &GetClassIdentifier() const
//    {
//      return GetStaticClassName( *this);
//    }
//
//    int Run() const
//    {
//
//    //////////////////////////////////
//    // construction and destruction //
//    //////////////////////////////////
//
//      // construct from properties
//      chemistry::FragmentRemoveBond fragment_remove_bond;
//
//      // clone function
//      util::ShPtr< chemistry::FragmentRemoveBond> fragment_remove_bond_clone( fragment_remove_bond.Clone());
//      BCL_ExampleCheck( fragment_remove_bond_clone.IsDefined(), true);
//
//      // construct with arguments
////      chemistry::FragmentRemoveBond fragment_remove_bonder
////      (
////        chemistry::CollectorValence(),
////        chemistry::PickAtomRandom()
////      );
//
//    /////////////////
//    // data access //
//    /////////////////
//
//    ///////////////
//    // operators //
//    ///////////////
//
//    ////////////////
//    // operations //
//    ////////////////
//
//      // create input stream for reading a base molecule
//      io::IFStream input;
//      BCL_ExampleMustOpenInputFile( input, AddExampleInputPathToFilename( e_Chemistry, "d6w.sdf"));
//
//      // read in ensemble
//      chemistry::FragmentEnsemble base_mol( input);
//      chemistry::FragmentComplete d6w( base_mol.GetMolecules().FirstElement());
//
//      // close stream
//      io::File::CloseClearFStream( input);
//
//      // Transform an atom in d6w
////      math::MutateResult< chemistry::FragmentComplete> new_mutate;
////      util::ShPtr< chemistry::FragmentComplete> new_molecule( new chemistry::FragmentComplete);
////      new_mutate = fragment_remove_bond( d6w);
////      BCL_Debug( new_mutate);
////      if( new_mutate.GetArgument().IsDefined())
////      {
////        new_molecule = new_mutate.GetArgument();
////        d6w = *new_molecule;
////      }
////
////        return 0;
////
////      // Write the Generated fragment to a file
////      const std::string fragment_remove_bond_out( AddExampleOutputPathToFilename( *new_molecule, "fragment_remove_bond.sdf"));
////      io::OFStream output;
////      BCL_ExampleMustOpenOutputFile( output, fragment_remove_bond_out);
////
////      // write out the constitution as we do not have coordinates
////      new_molecule->WriteMDL( output);
////      io::File::CloseClearFStream( output);
////      std::string correct_filename( fragment_remove_bond_out + ".correct");
//
//    //////////////////////
//    // input and output //
//    //////////////////////
//
//    //////////////////////
//    // helper functions //
//    //////////////////////
//
//      return 0;
//    } // Run
//
//    static const ExampleClass::EnumType s_Instance;
//
//  }; //end ExampleChemistryFragmentRemoveBond
//
//  const ExampleClass::EnumType ExampleChemistryFragmentRemoveBond::s_Instance
//  (
//    GetExamples().AddEnum( ExampleChemistryFragmentRemoveBond())
//  );
//
//} // namespace bcl
