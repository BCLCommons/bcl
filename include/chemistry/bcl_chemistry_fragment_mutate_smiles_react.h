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

#ifndef BCL_CHEMISTRY_FRAGMENT_MUTATE_SMILES_REACT_H_
#define BCL_CHEMISTRY_FRAGMENT_MUTATE_SMILES_REACT_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically
#include "descriptor/bcl_descriptor.fwd.hh"
#include "find/bcl_find.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_atom_conformational_interface.h"
#include "bcl_chemistry_collector_valence.h"
#include "bcl_chemistry_fragment_complete.h"
#include "bcl_chemistry_fragment_constitution_shared.h"
#include "bcl_chemistry_fragment_ensemble.h"
#include "bcl_chemistry_fragment_mutate_interface.h"
#include "descriptor/bcl_descriptor_base.h"
#include "find/bcl_find_pick_interface.h"
#include "math/bcl_math_mutate_interface.h"
#include "math/bcl_math_mutate_result.h"
#include "util/bcl_util_function_interface.h"
#include "util/bcl_util_sh_ptr.h"
#include "util/bcl_util_si_ptr_list.h"
// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class FragmentMutateSmilesReact
    //! @brief Used to read reaction SMILES/SMIRKS/SMARTS files to modify molecules during design
    //!
    //! @see @link example_chemistry_fragment_mutate_smiles_react.cpp @endlink
    //! @author brownbp1
    //! @date Dec 07, 2022
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API FragmentMutateSmilesReact :
      public FragmentMutateInterface
    {

    public:

      ////////////////////
      // helper classes //
      ////////////////////

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //!
      //! @class SmilesReactionComponent
      //! @brief A container object describing reagent components of SMIRKS reaction
      //!
      //! @author brownbp1
      //! @date Dec 7, 2022
      //!
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      struct SmilesReactionComponent {

        //! Reagent SMILES/SMARTS string
        std::string m_ReagentSmiles;

        //! Group ID for the reagent
        size_t m_ReagentID;

        //! Reaction ID
        std::string m_ReactionID;

        //! Position in the reaction that this reagent holds
        size_t m_ReactionPosition;

        // Constructors

        //! @brief Default constructor
        SmilesReactionComponent() :
          m_ReagentSmiles( ""),
          m_ReagentID( util::GetUndefinedSize_t()),
          m_ReactionID( ""),
          m_ReactionPosition( util::GetUndefinedSize_t())
        {
        }

        //! @brief Construct from all data
        SmilesReactionComponent
        (
          const std::string &SMILES,
          const size_t GROUP_ID,
          const std::string &RXN_ID,
          const size_t RXN_POS
        ) :
          m_ReagentSmiles( SMILES),
          m_ReagentID( GROUP_ID),
          m_ReactionID( RXN_ID),
          m_ReactionPosition( RXN_POS)
        {
        }

        // Operators
        // @brief less-than operator
        // @return true if the lhs object is less than the rhs object; false otherwise
        bool operator<( SmilesReactionComponent const &RHS)
        {
          if( m_ReagentID < RHS.m_ReagentID )
          {
            return true;
          }
          return false;
        }
      };

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //!
      //! @class SmirksReactor
      //! @brief A poorly written, incomplete parser for SMIRKS reaction files
      //!
      //! @author brownbp1
      //! @date Dec 8, 2022
      //!
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      struct SmirksReactor {

        //! @brief full reaction data string
        //! @details consists of 8 space-separated columns:
        //! reaction_id number_rxn_positions smirks_reaction smirks_product smirks_reagent_1 smirks_reagent_2 smirks_reagent_3 smirks_reagent_4
        std::string m_ReactionData;

        //! @brief reaction ID
        std::string m_SmirksReactionID;

        //! @brief reaction
        std::string m_SmirksReaction;

        //! @brief number of reagents (also the number of reaction positions)
        size_t m_NumberReagents;

        //! @brief reagents
        storage::Vector< std::string> m_SmirksReagents;


        //! @brief Return the BCL bond type from a SMIRKS bond type
        //! @details SMIRKS bond types are simply encoded for single, double,
        //! and triple bonds, as well as an ambiguous "any" bond type and a generic
        //! "ring" bond type. The BCL bond types are higher resolution, and there
        //! may be multiple that fit under the SMIRKS bond type umbrella; however,
        //! this function returns the most simple and generic version of the bond
        //! type under the hope that if molecule perturbations or conversions
        //! are occurring that the molecule is standardized after-the-fact
        //! @return the BCL configurational bond type
        ConfigurationalBondType GetBondTypeFromSmirks( std::string &SMIRKS_BOND_TYPE)
        {
          if( SMIRKS_BOND_TYPE == "-")
          {
            return GetConfigurationalBondTypes().e_NonConjugatedSingleBond;
          }
          else if( SMIRKS_BOND_TYPE == "=")
          {
            return GetConfigurationalBondTypes().e_ConjugatedDoubleBond;
          }
          else if( SMIRKS_BOND_TYPE == "#")
          {
            return GetConfigurationalBondTypes().e_ConjugatedTripleBond;
          }
          else if( SMIRKS_BOND_TYPE == "~")
          {
            return GetConfigurationalBondTypes().e_Undefined; // TODO reactions involving these will require more sophisticated parsing
          }
          else if( SMIRKS_BOND_TYPE == "@")
          {
            return GetConfigurationalBondTypes().e_ConjugatedBondInRing; // consider also e_NonConjugatedSingleBondInRing
          }
          else
          {
            return GetConfigurationalBondTypes().e_Undefined; // TODO conflicting with ~
          }
        }


        // TODO this really needs to be private when this struct is refactored into its own class
        //! @brief Get dummy/reactant atom from reagent string
        storage::Vector< ElementType> ParseDummyAtom( const std::string &REAGENT)
        {
          // split apart the components of the reagent
          BCL_Debug( REAGENT);
          storage::Vector< std::string> reagent_components( util::SplitString( REAGENT, "[]")); // TODO verify the split
          BCL_Debug( reagent_components);

          // iterate over components and collect the valid element types that pop out
          storage::Vector< ElementType> elements;
          for
          (
              auto itr( reagent_components.Begin()), itr_end( reagent_components.End());
              itr != itr_end;
              ++itr
          )
          {
            // skip bond type strings
            if
            (
                *itr == "-" ||
                *itr == "=" ||
                *itr == "#" ||
                *itr == "~" ||
                *itr == "@"
            )
            {
              continue;
            }

            ElementType element( GetElementTypes().ElementTypeLookup( *itr));
            if( element != GetElementTypes().e_Undefined)
            {
              elements.PushBack( element);
            }
            else
            {
              BCL_MessageStd( "Skipping undefined element!");
            }
          }
          BCL_Debug( elements);
          return elements;
        }

        // TODO this really needs to be private when this struct is refactored into its own class
        //! @brief Get dummy/reactant atom from reagent string
        storage::Map< ElementType, ConfigurationalBondType> ParseDummyAtomBonds( const std::string &REAGENT)
        {
          // initialize output
          storage::Map< ElementType, ConfigurationalBondType> mapped_bonds;

          // get our dummy atoms
          storage::Vector< ElementType> dummy_elements( ParseDummyAtom( REAGENT ) );
          if( !dummy_elements.GetSize())
          {
            return storage::Map< ElementType, ConfigurationalBondType>();
          }

          // for each dummy atom check its bond to a neighbor atom
          for( size_t e_i( 0), e_sz( dummy_elements.GetSize()); e_i < e_sz; ++e_i)
          {
            const ElementType &element_type( dummy_elements( e_i));
            BCL_Debug( element_type);
            std::string smirks_ele( "[" +  element_type->GetChemicalSymbol() + "]");
            BCL_Debug( smirks_ele);

            // find where this dummy atom occurs in the original reagent string
            // assumptions: (1) getting first occurrence is sufficient;
            // (2) if there is a second occurrence then the bond type would be the same as it indicates
            // a cyclical closure within the reagent
            size_t dummy_atom_npos( REAGENT.find( smirks_ele));
            BCL_Debug( dummy_atom_npos);
            size_t reagent_str_size( REAGENT.size());
            BCL_Debug( reagent_str_size);

            // left-edge: only check +1 npos
            ConfigurationalBondType bond_type;
            if( dummy_atom_npos == size_t( 0))
            {
              std::string bond( REAGENT.substr( dummy_atom_npos + 3, 1));
              BCL_Debug( bond);
              bond_type = this->GetBondTypeFromSmirks( bond);
            }
            // right-edge; only check -1 npos
            else if( dummy_atom_npos == reagent_str_size - 1)
            {
              std::string bond( REAGENT.substr( dummy_atom_npos - 1, 1));
              BCL_Debug( bond);
              bond_type = this->GetBondTypeFromSmirks( bond);
            }
            // otherwise check both -1 and +1 npos
            else
            {
              std::string bond( REAGENT.substr( dummy_atom_npos - 1, 1));
              BCL_Debug( bond);
              bond_type = this->GetBondTypeFromSmirks( bond);

              if( bond_type == GetConfigurationalBondTypes().e_Undefined) // TODO generally, all of this is incompatible with '~' bond types
              {
                bond = REAGENT.substr( dummy_atom_npos + 3, 1);
                bond_type = this->GetBondTypeFromSmirks( bond);
              }
            }
            mapped_bonds.Insert( storage::Pair< ElementType, ConfigurationalBondType>( element_type, bond_type));
          }
          BCL_Debug( mapped_bonds);
          return mapped_bonds;
        }


        //! @brief Parses reaction data
        void ParseReactionData( const std::string &DATA)
        {
          // first split the line into the 8 columns
          BCL_Debug( DATA);
          storage::Vector< std::string> split_data( util::SplitString(DATA, " \t"));
          BCL_Debug( split_data);

          // convenience
          m_SmirksReactionID = split_data( 0);
          m_NumberReagents = util::ConvertStringToNumericalValue< size_t>( split_data( 1) );
          m_SmirksReaction = split_data( 2);

          // now split the reaction taking the first n_reagents positions as reagents
          BCL_Debug( m_SmirksReaction);
          storage::Vector< std::string> smirks_reagents_products( util::SplitString(m_SmirksReaction, ".'>>'")); // TODO this seems hairy
          BCL_Debug( smirks_reagents_products);
          storage::Vector< std::string> smirks_reagents;
          for( size_t i( 0); i < m_NumberReagents; ++i)
          {
            smirks_reagents.PushBack( smirks_reagents_products( i));
          }
          m_SmirksReagents = smirks_reagents;
//          const std::string &smirks_products( smirks_reagents_products( m_NumberReagents));

          // make sure that we have the correct number of reagents and products
          BCL_Assert
          (
            smirks_reagents_products.GetSize() == smirks_reagents.GetSize() + 1, // only ever expect 1 product
            "[ERROR] SmirksReactor::ParseReactionData reagents and products incorrectly parsed!"
          );

        }

        //! @brief Parses reaction data
        void ParseReactionData()
        {
          ParseReactionData( m_ReactionData);
        }


        //! @brief default constructor
        SmirksReactor() :
          m_ReactionData( ""),
          m_SmirksReaction( ""),
          m_NumberReagents( 0),
          m_SmirksReagents( storage::Vector< std::string>())
        {
        }

        //! @brief data constructor
        SmirksReactor( const std::string &DATA) :
          m_ReactionData( DATA)
        {
          ParseReactionData();
        }

      };

    /////////////
    // friends //
    /////////////

    private:

      //! @brief allowed reactions
      std::string m_ReactionFilename;
      std::string m_ReactionFileContents;
      storage::Vector< std::string> m_ReactionIDs;
      std::map< std::string, std::string> m_Reactions;
      bool m_InitializedReactions;

      //! @brief allowed reagents with which to perform reactions
      std::string m_ReagentsFilename;
      std::string m_ReagentsFileContents;
      bool m_InitializedReagents;

      //! @brief associated reagents and reactions;
      //! @details the key is a pair consisting of the reaction ID string and the reaction position;
      //! the value is a vector of all the SMILES reaction data containers matching the key
      std::map< std::pair< std::string, size_t>, std::vector< SmilesReactionComponent> > m_AssociatedReactions;

      //! @brief allowed reaction positions of starting molecule
      std::string m_AllowedRxnPositions;
      storage::Vector< size_t> m_AllowedRxnPositionIndices;

    //////////
    // data //
    //////////

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      FragmentMutateSmilesReact();

      //! @brief construct with a pool of external fragments for fragment grow
      //! @param FRAGMENT_POOL external fragments to add to base fragment
      FragmentMutateSmilesReact
      (
        const bool &CORINA_CONFS
      );

      //! @brief druglikeness constructor
      //! @param FRAGMENT_POOL external fragments to add to base fragment
      //! @param DRUG_LIKENESS_TYPE type of druglikeness filter to apply during clean
      FragmentMutateSmilesReact
      (
        const std::string &DRUG_LIKENESS_TYPE,
        const bool &CORINA_CONFS
      );

      //! @brief local mutate constructor
      //! @param FRAGMENT_POOL external fragments to add to base fragment
      //! @param DRUG_LIKENESS_TYPE type of druglikeness filter to apply during clean
      //! @param SCAFFOLD_FRAGMENT fragment to which the new mutated molecule will be aligned based on substructure
      //! @param MUTABLE_FRAGMENTS non-mutable component of the current molecule
      //! @param MUTABLE_ATOM_INDICES indices of atoms that can be mutated
      FragmentMutateSmilesReact
      (
        const std::string &DRUG_LIKENESS_TYPE,
        const FragmentComplete &SCAFFOLD_FRAGMENT,
        const FragmentEnsemble &MUTABLE_FRAGMENTS,
        const storage::Vector< size_t> &MUTABLE_ATOM_INDICES,
        const bool &CORINA_CONFS
      );

      //! @brief local mutate pose-sensitive constructor
      //! @param FRAGMENT_POOL external fragments to add to base fragment
      //! @param DRUG_LIKENESS_TYPE type of druglikeness filter to apply during clean
      //! @param SCAFFOLD_FRAGMENT fragment to which the new mutated molecule will be aligned based on substructure
      //! @param MUTABLE_FRAGMENTS non-mutable component of the current molecule
      //! @param MUTABLE_ATOM_INDICES indices of atoms that can be mutated
      //! @param MDL property label containing path to protein binding pocket PDB file
      //! @param PROPERTY_SCORER property that will be used to score interactions with protein pocket
      //! @param RESOLVE_CLASHES if true, resolve clashes with specified protein pocket after mutatation
      //! @param BFACTORS vector of values indicating per-residue flexibility (higher values are more flexible)
      FragmentMutateSmilesReact
      (
        const std::string &DRUG_LIKENESS_TYPE,
        const FragmentComplete &SCAFFOLD_FRAGMENT,
        const FragmentEnsemble &MUTABLE_FRAGMENTS,
        const storage::Vector< size_t> &MUTABLE_ATOM_INDICES,
        const std::string &MDL,
        const descriptor::CheminfoProperty &PROPERTY_SCORER,
        const bool &RESOLVE_CLASHES,
        const storage::Vector< float> &BFACTORS,
        const bool &CORINA_CONFS
      );

      //! @brief local clash resolver constructor
      //! @param FRAGMENT_POOL external fragments to add to base fragment
      //! @param DRUG_LIKENESS_TYPE type of druglikeness filter to apply during clean
      //! @param SCAFFOLD_FRAGMENT fragment to which the new mutated molecule will be aligned based on substructure
      //! @param MUTABLE_FRAGMENTS non-mutable component of the current molecule
      //! @param MUTABLE_ATOM_INDICES indices of atoms that can be mutated
      //! @param MDL property label containing path to protein binding pocket PDB file
      //! @param RESOLVE_CLASHES if true, resolve clashes with specified protein pocket after mutatation
      //! @param BFACTORS vector of values indicating per-residue flexibility (higher values are more flexible)
      FragmentMutateSmilesReact
      (
        const std::string &DRUG_LIKENESS_TYPE,
        const FragmentComplete &SCAFFOLD_FRAGMENT,
        const FragmentEnsemble &MUTABLE_FRAGMENTS,
        const storage::Vector< size_t> &MUTABLE_ATOM_INDICES,
        const std::string &MDL,
        const bool &RESOLVE_CLASHES,
        const storage::Vector< float> &BFACTORS,
        const bool &CORINA_CONFS
      );

      //! @brief clone constructor
      FragmentMutateSmilesReact *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief returns the name used for this class in an object data label
      //! @return the name used for this class in an object data label
      const std::string &GetAlias() const;

    ///////////////
    // operators //
    ///////////////

      //! @brief virtual operator taking an SmallMolecule and returning a grown SmallMolecule
      //! @param FRAGMENT small molecule of interest
      //! @return Constitution after the mutate
      math::MutateResult< FragmentComplete> operator()( const FragmentComplete &FRAGMENT) const;

    ////////////////
    // operations //
    ////////////////

      //! @brief set medchem fragment library from filename
      void SetFragmentLibraryFromFilename( const std::string &FRAGMENTS_FILENAME);

      //! @brief get a random reagent matching the reaction id and position
      //! @return a reagent molecule
      FragmentComplete GetRandomReagent( const std::string &REACTION_ID, const size_t REACTION_POS) const;

      //! @brief get a random reaction
      //! @return a reaction ID string
      std::string GetRandomReactionID() const;

      //! @brief get a random reaction
      //! @return a reaction ID string
      std::string GetReactionFromID( const std::string &RXN_ID) const;

      //! @brief get random allowed reaction position from user-specified options
      size_t GetRandomReactionPosition( const storage::Vector< size_t> &RXN_POSITIONS) const;

    private:

      //! @brief combine two fragments through their pseudoreaction scheme
      //! @details attach two fragments using their mutually matched dummy atom element types
      //! and mapped attachment atom indices. After attaching the two fragments, it is
      //! possible that the new product molecule contains additional mutually matched
      //! dummy elements within a single fragment. Therefore, after combining the two
      //! fragments, this function will perform intramolecular pseudoreactions until
      //! there are no more mutually matched dummy atoms
      //! @return a product molecule with an associated map between un-reacted dummy atom
      //! element types and the attachment indices using the new product molecule indexing
      storage::Pair< FragmentComplete, storage::Map< ElementType, size_t>> ReactFragments
      (
        const storage::Pair< FragmentComplete, storage::Map< ElementType, size_t>> &REAGENT_A,
        const storage::Pair< FragmentComplete, storage::Map< ElementType, size_t>> &REAGENT_B,
        const storage::Map< ElementType, ConfigurationalBondType> &BONDS
      ) const;

      //! @brief recursively call ReactFragments to combine reagents into a single product
      //! @return a product molecule with an associated map between un-reacted dummy atom
      //! element types and the attachment indices using the new product molecule indexing
      storage::Pair< FragmentComplete, storage::Map< ElementType, size_t>> ReactFragments
      (
        const storage::Vector< storage::Pair< FragmentComplete, storage::Map< ElementType, size_t>>> &REAGENTS,
        const storage::Map< ElementType, ConfigurationalBondType> &BONDS
      ) const;

      //! @brief remove the dummy element from the input molecule
      //! @return the new molecule with all dummy elements removed, as well
      //! as a map between the dummy element type and the atom that used to connect
      //! to that dummy element
      storage::Pair< FragmentComplete, storage::Map< ElementType, size_t>> RemoveDummyElement
      (
        const FragmentComplete &MOLECULE,
        const storage::Vector< ElementType> &ELEMENTS
      ) const;

    //////////////////////
    // helper functions //
    //////////////////////

    protected:

      //! @brief Set the data members corresponding to reagents files
      //! @return true if initialization successful; false otherwise
      bool InitializeReagents();

      //! @brief Set the data members corresponding to reactions files
      //! @return true if initialization successful; false otherwise
      bool InitializeReactions();

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

      //! @brief Set the members of this property from the given LABEL
      //! @param LABEL the label to parse
      //! @param ERROR_STREAM the stream to write errors to
      bool ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERROR_STREAM);

    //////////////////////
    // input and output //
    //////////////////////

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM output stream
      //! @param INDENT number of indentations
      //! @return ostream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; // class FragmentMutateSmilesReact

  } // namespace chemistry
} // namespace bcl

#endif //BCL_CHEMISTRY_FRAGMENT_MUTATE_SMILES_REACT_H_
