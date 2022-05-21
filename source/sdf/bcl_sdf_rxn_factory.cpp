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

// initialize the static initialization fiasco finder, if macro ENABLE_FIASCO_FINDER is defined
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "sdf/bcl_sdf_rxn_factory.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_bond_isometry_handler.h"
#include "chemistry/bcl_chemistry_fragment_complete.h"
#include "chemistry/bcl_chemistry_reaction_complete.h"
#include "sdf/bcl_sdf_fragment_factory.h"

namespace bcl
{
  namespace sdf
  {

    //! @brief construct a ReactionComplete from a rxn file
    //! @param HANDLER the handler to use
    //! @param H_PREF hydrogen handling preference
    chemistry::ReactionComplete RXNFactory::MakeReactionComplete( const RXNHandler &HANDLER)
    {
      // The MDL handlers for products and reactants
      //const storage::Vector< MdlHandler> &reactant_handlers( HANDLER.GetReactantMdlHandlers());
      //const storage::Vector< MdlHandler> &product_handlers( HANDLER.GetProductMdlHandlers());
      const storage::Vector< MolfileHandler> &reactant_handlers( HANDLER.GetReactantHandlers());
      const storage::Vector< MolfileHandler> &product_handlers( HANDLER.GetProductHandlers());
      storage::Vector< storage::Vector< AtomInfo> > reactant_atom_info;
      storage::Vector< storage::Vector< BondInfo> > reactant_bond_info;

      storage::Vector< storage::Vector< AtomInfo> > product_atom_info;
      storage::Vector< storage::Vector< BondInfo> > product_bond_info;

      //storage::Vector< chemistry::FragmentComplete> reactants;
      //storage::Vector< chemistry::FragmentComplete> products;

      // make the reactants
      storage::Vector< chemistry::ReactionStructure> reactant_structs;
      for( size_t r( 0); r < reactant_handlers.GetSize(); ++r)
      {
        reactant_structs.PushBack( MakeReactionStructure( reactant_handlers( r)));
      }

      // Make the products
      storage::Vector< chemistry::ReactionStructure> product_structs;
      for( size_t p( 0); p < product_handlers.GetSize(); ++p)
      {
        product_structs.PushBack( MakeReactionStructure( product_handlers( p)));
      }

      // Make the reaction complete
      return chemistry::ReactionComplete
      (
        reactant_structs,
        product_structs,
        HANDLER.GetReactiveAtomsInReactants(),
        HANDLER.GetReactiveAtomsInProducts(),
        HANDLER.GetDescription()
      );

    }

    //! @brief construct a ReactionStructure from a molfile
    chemistry::ReactionStructure RXNFactory::MakeReactionStructure( const CTabHandler &CTAB)
    {
      chemistry::ReactionStructure rxn_struct;
      if( !CTAB.IsValid())
      {
        return rxn_struct;
      }

      const storage::Map< std::string, storage::Vector< std::string> > &mdl_props( CTAB.GetCachedProperties());

      chemistry::AtomVector< chemistry::AtomComplete> atom_v( CTAB.GetAtomInfo(), CTAB.GetBondInfo());
      chemistry::FragmentComplete new_frag( atom_v, std::string());
      
      std::vector< chemistry::ReactionStructure::Aromaticity> allowed_aromaticities( atom_v.GetSize(), chemistry::ReactionStructure::e_AliphaticOrAromatic);

      const std::string arom_prop( MdlProperty( MdlProperty::e_BclAtomAromaticity).GetLabel());
      if( mdl_props.Has( arom_prop))
      {
        const storage::Vector< std::string> &arom_str( mdl_props.GetValue( arom_prop));
        std::vector< size_t> arom_atoms( atom_v.GetSize(), size_t( 0));

        // store the aromaticity of the molecule
        for( size_t i( 0); i < arom_str.GetSize(); i += 2)
        {
          size_t atom_no( 0);
          if( !util::TryConvertFromString( atom_no, arom_str( i), util::GetLogger()))
          {
            BCL_MessageStd( "Could not parse aromaticity information from reaction");
            allowed_aromaticities = std::vector< chemistry::ReactionStructure::Aromaticity>( atom_v.GetSize(), chemistry::ReactionStructure::e_AliphaticOrAromatic);
            break;
          }
          if( atom_no == 0 || atom_no > atom_v.GetSize())
          {
            BCL_MessageStd( "Could not parse aromaticity information from reaction");
            if( !atom_no)
            {
              BCL_MessageStd( "Atom number was 0, minimum value is 1 as per the SDF numbering standard");
            }
            allowed_aromaticities = std::vector< chemistry::ReactionStructure::Aromaticity>( atom_v.GetSize(), chemistry::ReactionStructure::e_AliphaticOrAromatic);
            break;
          }
          if( arom_str( i + 1) == "ARO")
          {
            allowed_aromaticities[ atom_no - 1] = chemistry::ReactionStructure::e_Aromatic;
          }
          else if( arom_str( i + 1) == "ALI")
          {
            allowed_aromaticities[ atom_no - 1] = chemistry::ReactionStructure::e_Aliphatic;
          }
        }
      }
      
      rxn_struct = chemistry::ReactionStructure( new_frag, allowed_aromaticities);

      return rxn_struct;
    }

  } // namespace sdf
} // namespace bcl

