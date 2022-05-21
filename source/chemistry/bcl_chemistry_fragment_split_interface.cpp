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
#include "chemistry/bcl_chemistry_fragment_split_interface.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

  ////////////////
  // operations //
  ////////////////

    //! @brief splits the molecule according to GetComponentVertices
    //! @param CONFORMATION the molecule to split
    FragmentEnsemble FragmentSplitInterface::operator()( const ConformationInterface &CONFORMATION) const
    {
      // create a graph
      ConformationGraphConverter::t_AtomGraph mol_graph
      (
        ConformationGraphConverter::CreateGraphWithAtoms( CONFORMATION)
      );
      // get connected componenets of graph
      storage::List< storage::Vector< size_t> > component_vertices( GetComponentVertices( CONFORMATION, mol_graph));

      // return an ensemble of fragments derived from componenets of graph
      return ConvertComponentsIntoEnsemble( CONFORMATION, component_vertices, mol_graph);
    }

    //! @brief helper function that can be used to convert the output of GetComponentVertices into an ensemble
    //! @param MOLECULE the molecule to split
    //! @param COMPONENTS list of components to create
    //! @param MOLECULE_GRAPH graph of the molecule with atoms
    //! @return a fragment ensemble
    FragmentEnsemble FragmentSplitInterface::ConvertComponentsIntoEnsemble
    (
      const ConformationInterface &MOLECULE,
      const storage::List< storage::Vector< size_t> > &COMPONENTS,
      const ConformationGraphConverter::t_AtomGraph &MOLECULE_GRAPH,
      const bool &INVERT_OUTPUT
    ) const
    {
      // ensemble to store components of molecule of interest
      FragmentEnsemble ensemble_of_components;

      // atom indices of the current molecule
      storage::Vector< size_t> atom_indices;
      for
      (
          auto atom_itr( MOLECULE.GetAtomsIterator().Begin()), atom_itr_end( MOLECULE.GetAtomsIterator().End());
          atom_itr != atom_itr_end;
          ++atom_itr
      )
      {
        atom_indices.PushBack( MOLECULE.GetAtomIndex( *atom_itr));
      }

      // add all split components to the ensemble
      for
      (
        storage::List< storage::Vector< size_t> >::const_iterator
          itr( COMPONENTS.Begin()), itr_end( COMPONENTS.End());
        itr != itr_end;
        ++itr
      )
      {
        // if size of component is greater than minimum size then continue
        if( itr->GetSize() > GetMinSize())
        {
          // add component to output ensemble
          if( !INVERT_OUTPUT)
          {
            FragmentComplete new_fragment
            (
              ConformationGraphConverter::CreateAtomsFromGraph( MOLECULE_GRAPH.GetSubgraph( *itr)),
              GetAlias() + " isolated from " + MOLECULE.GetName()
            );
            ensemble_of_components.PushBack( new_fragment);
          }
          // add everything except component to output ensemble
          else
          {
            storage::Vector< size_t> keep_atom_indices( atom_indices);
            for
            (
                auto comp_index_itr( itr->Begin()), comp_index_itr_end( itr->End());
                comp_index_itr != comp_index_itr_end;
                ++comp_index_itr
            )
            {
              auto atom_pos_itr( std::find( keep_atom_indices.Begin(), keep_atom_indices.End(), *comp_index_itr));
              if( atom_pos_itr != keep_atom_indices.End())
              {
                size_t pos( std::distance( keep_atom_indices.Begin(), atom_pos_itr));
                keep_atom_indices.RemoveElements( pos, 1);
              }
            }
            std::unique( keep_atom_indices.Begin(), keep_atom_indices.End());
            FragmentComplete new_fragment
            (
              ConformationGraphConverter::CreateAtomsFromGraph( MOLECULE_GRAPH.GetSubgraph( keep_atom_indices)),
              GetAlias() + " isolated from " + MOLECULE.GetName()
            );
            ensemble_of_components.PushBack( new_fragment);
          }
        }
      }
      return ensemble_of_components;
    }

  } // namespace chemistry
} // namespace bcl
