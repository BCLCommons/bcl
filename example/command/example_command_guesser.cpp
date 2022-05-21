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
#include "command/bcl_command_guesser.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_command_guesser.cpp
  //!
  //! @author mendenjl
  //! @date Nov 17, 2012
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleCommandGuesser :
    public ExampleInterface
  {
  public:

    ExampleCommandGuesser *Clone() const
    {
      return new ExampleCommandGuesser( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! returns class name
    //! the class name as const ref std::string
    const std::string &GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    int Run() const
    {

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // create a vector with a few different hard-to-remember options
      storage::Vector< std::string> possible_values
      (
        storage::Vector< std::string>::Create
        (
          "this is my favorite",
          "OptionToMispell",
          "Generate",
          "Generator",
          "Generation",
          "Generates",
          "ALPHA   BETA   GAMMA",
          "ChemistryMolecule",
          "ChemistryMoleculeProperties",
          "ChemistryMoleculePropertiesBondTypeCount"
        )
      );

      // add some other flags of interest
      possible_values.Append
      (
        storage::Vector< std::string>::Create
        (
          "AA",
          "Res",
          "aa_ss_pred_list"
        )
      );

      // ask the default guesser to try guess what we meant
      BCL_ExampleCheck
      (
        command::Guesser::GetDefaultGuesser().Guess( "yep this is my fave", possible_values).First(),
        command::Guesser::e_SomeWordsMatch
      );

      BCL_ExampleCheck
      (
        command::Guesser::GetDefaultGuesser().Guess( "yep this is my fave", possible_values).Second(),
        storage::Vector< std::string>( size_t( 1), possible_values( 0))
      );

      BCL_ExampleCheck
      (
        command::Guesser::GetDefaultGuesser().Guess( "favorite is my this", possible_values).First(),
        command::Guesser::e_ReorderedWords
      );

      BCL_ExampleCheck
      (
        command::Guesser::GetDefaultGuesser().Guess( "favorite is my this", possible_values).Second(),
        storage::Vector< std::string>( size_t( 1), possible_values( 0))
      );

      BCL_ExampleCheck
      (
        command::Guesser::GetDefaultGuesser().Guess( "this i", possible_values).First(),
        command::Guesser::e_FirstLetters
      );

      BCL_ExampleCheck
      (
        command::Guesser::GetDefaultGuesser().Guess( "favorites is my this", possible_values).First(),
        command::Guesser::e_ReorderedStems
      );

      BCL_ExampleCheck
      (
        command::Guesser::GetDefaultGuesser().Guess( "this i", possible_values).Second(),
        storage::Vector< std::string>( size_t( 1), possible_values( 0))
      );

      BCL_ExampleCheck
      (
        command::Guesser::GetDefaultGuesser().Guess( "alpha beta gamma", possible_values).First(),
        command::Guesser::e_CaseOrSpace
      );

      BCL_ExampleCheck
      (
        command::Guesser::GetDefaultGuesser().Guess( "alpha beta gamma", possible_values).Second(),
        storage::Vector< std::string>( size_t( 1), possible_values( 6))
      );

      BCL_ExampleCheck
      (
        command::Guesser::GetDefaultGuesser().Guess( "OptionToMispells", possible_values).First(),
        command::Guesser::e_Suffix
      );

      BCL_ExampleCheck
      (
        command::Guesser::GetDefaultGuesser().Guess( "OptionToMispells", possible_values).Second(),
        storage::Vector< std::string>( size_t( 1), possible_values( 1))
      );

      BCL_ExampleCheck
      (
        command::Guesser::GetDefaultGuesser().Guess( "Chem", possible_values).First(),
        command::Guesser::e_FirstLetters
      );
      BCL_ExampleCheck
      (
        command::Guesser::GetDefaultGuesser().Guess( "Chem", possible_values).Second().GetSize(),
        3
      );

      BCL_ExampleCheck
      (
        command::Guesser::GetDefaultGuesser().Guess( "ChemMol", possible_values).First(),
        command::Guesser::e_StrongAbbreviation
      );

      BCL_ExampleCheck
      (
        command::Guesser::GetDefaultGuesser().Guess( "ChemMol", possible_values).Second(),
        storage::Vector< std::string>::Create
        (
          "ChemistryMolecule",
          "ChemistryMoleculeProperties",
          "ChemistryMoleculePropertiesBondTypeCount"
        )
      );

      BCL_ExampleCheck
      (
        command::Guesser::GetDefaultGuesser().Guess( "CMP", possible_values).First(),
        command::Guesser::e_StrongAbbreviation
      );

      BCL_ExampleCheck
      (
        command::Guesser::GetDefaultGuesser().Guess( "CMP", possible_values).Second(),
        storage::Vector< std::string>::Create
        (
          "ChemistryMoleculeProperties",
          "ChemistryMoleculePropertiesBondTypeCount"
        )
      );

      BCL_ExampleCheck
      (
        command::Guesser::GetDefaultGuesser().Guess( "Cmp", possible_values).First(),
        command::Guesser::e_WeakAbbreviation
      );

      // try something that with only some of the capital letters in the full abreviation
      BCL_ExampleCheck
      (
        command::Guesser::GetDefaultGuesser().Guess( "CMPBT", possible_values).First(),
        command::Guesser::e_StrongAbbreviation
      );
      BCL_ExampleCheck
      (
        command::Guesser::GetDefaultGuesser().Guess( "CMPBT", possible_values).Second(),
        storage::Vector< std::string>( size_t( 1), "ChemistryMoleculePropertiesBondTypeCount")
      );

      // try a longer abbreviation
      BCL_ExampleCheck
      (
        command::Guesser::GetDefaultGuesser().Guess( "CMPBTC", possible_values).First(),
        command::Guesser::e_StrongAbbreviation
      );

      // try spelling words out
      BCL_ExampleCheck
      (
        command::Guesser::GetDefaultGuesser().Guess( "AminoAcid", possible_values).First(),
        command::Guesser::e_Explicit
      );
      BCL_ExampleCheck
      (
        command::Guesser::GetDefaultGuesser().Guess( "AminoAcid", possible_values).Second(),
        storage::Vector< std::string>( size_t( 1), "AA")
      );

      // should work even if the full word is partially abbreviated
      BCL_ExampleCheck
      (
        command::Guesser::GetDefaultGuesser().Guess( "Residue", possible_values).Second(),
        storage::Vector< std::string>( size_t( 1), "Res")
      );

      // try it with flags
      BCL_ExampleCheck
      (
        command::Guesser::GetDefaultGuesser().Guess( "AASSPredList", possible_values).Second(),
        storage::Vector< std::string>( size_t( 1), "aa_ss_pred_list")
      );

      // try something else that is pretty ridiculous
      BCL_ExampleCheck
      (
        command::Guesser::GetDefaultGuesser().Guess( "amino_acid_secondary_structure_prediction_list", possible_values).Second(),
        storage::Vector< std::string>( size_t( 1), "aa_ss_pred_list")
      );

      // try using punctuation
      BCL_ExampleCheck
      (
        command::Guesser::GetDefaultGuesser().Guess
        (
          "const(0)",
          storage::Vector< std::string>( size_t( 1), "Constant")
        ).Second(),
        storage::Vector< std::string>( size_t( 1), "Constant")
      );

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleCommandGuesser

  const ExampleClass::EnumType ExampleCommandGuesser::s_Instance
  (
    GetExamples().AddEnum( ExampleCommandGuesser())
  );

} // namespace bcl

