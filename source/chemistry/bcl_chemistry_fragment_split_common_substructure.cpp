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
#include "chemistry/bcl_chemistry_fragment_split_common_substructure.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_conformation_graph_converter.h"
#include "chemistry/bcl_chemistry_fragment_complete.h"
#include "chemistry/bcl_chemistry_fragment_split_isolate.h"
#include "graph/bcl_graph_common_subgraph_isomorphism.h"
#include "graph/bcl_graph_connectivity.h"
#include "graph/bcl_graph_const_graph.h"
#include "io/bcl_io_file.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

    // add the interface to the set of known implementations
    const util::SiPtr< const util::ObjectInterface> FragmentSplitCommonSubstructure::s_Instance
    (
      util::Enumerated< FragmentSplitInterface>::AddInstance( new FragmentSplitCommonSubstructure)
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor, sets default steps to 4
    FragmentSplitCommonSubstructure::FragmentSplitCommonSubstructure
    (
      const std::string &FILENAME,
      const ConformationGraphConverter::AtomComparisonType ATOM_COMPARISON,
      const ConfigurationalBondTypeData::Data BOND_COMPARISON,
      const bool &COMPLEMENT
    ) :
      m_Converter( ATOM_COMPARISON, BOND_COMPARISON, false),
      m_File( FILENAME),
      m_AtomComparison( ATOM_COMPARISON),
      m_BondComparison( BOND_COMPARISON),
      m_Complement( COMPLEMENT),
      m_FileWasRead( false)
    {
    }

    //! virtual copy constructor
    FragmentSplitCommonSubstructure *FragmentSplitCommonSubstructure::Clone() const
    {
      return new FragmentSplitCommonSubstructure( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! the class name as const ref std::string
    const std::string &FragmentSplitCommonSubstructure::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &FragmentSplitCommonSubstructure::GetAlias() const
    {
      static const std::string s_name( "LargestCommonSubstructure");
      return s_name;
    }

    //! @brief Get a description for what this class does (used when writing help)
    //! @return a description for what this class does (used when writing help)
    const std::string &FragmentSplitCommonSubstructure::GetClassDescription() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gets the minimum size of fragments
    //! @return the minimum size of fragments
    const size_t FragmentSplitCommonSubstructure::GetMinSize() const
    {
      return 0;
    }

  /////////////////
  //  operations //
  /////////////////

    //! @brief returns an ensemble of fragments of a molecule
    //! @param CONFORMATION molecule of interest
    //! @return an ensemble of common substructures relative to those in a file
    //! TODO: Implement this
    FragmentEnsemble FragmentSplitCommonSubstructure::operator()( const ConformationInterface &CONFORMATION) const
    {

      // read in molecules (if necessary)
      ReadFile();

      graph::ConstGraph< size_t, size_t> mol_graph( m_Converter( CONFORMATION));
      util::OwnPtr< graph::ConstGraph< size_t, size_t> > mol_graph_ptr( &mol_graph, false);

      FragmentEnsemble common_structs;

      graph::CommonSubgraphIsomorphism< size_t, size_t> common_subgraph_iso;
      common_subgraph_iso.SetGraphA( mol_graph_ptr);
      for( size_t mol_no( 0); mol_no < m_MoleculeGraphs.GetSize(); ++mol_no)
      {
        // make a new version of the starting molecule
        FragmentComplete start_mol( CONFORMATION);

        util::OwnPtr< graph::ConstGraph< size_t, size_t> > mol_comp_graph_ptr( &m_MoleculeGraphs( mol_no), false);
        common_subgraph_iso.SetGraphB( mol_comp_graph_ptr);

        // get the actual isomorphism, using estimated upper bounds on its size
        common_subgraph_iso.FindIsomorphism( common_subgraph_iso.EstimateUpperBounds(), 1); // get all isomorphisms

        graph::Subgraph< size_t, size_t> subgraph
        (
          common_subgraph_iso.GetSubgraphIsomorphismsOfGraphA().FirstElement() // just returning the first iso? hmm...
        );

        // convert subgraph into small molecule
        AtomVector< AtomComplete> new_atom_vec( start_mol.GetAtomVector());
        if( m_Complement)
        {
          // BCL_Debug( subgraph.GetComplement().GetVertexIndices());
          new_atom_vec.Reorder( subgraph.GetComplement().GetVertexIndices());
        }
        else
        {
          // BCL_Debug( subgraph.GetVertexIndices());
          new_atom_vec.Reorder( subgraph.GetVertexIndices());
        }
        FragmentComplete split_mol( new_atom_vec, "");

        // if we have disconnected elements then we want those to be separate entries in the final SDF
        FragmentSplitIsolate isolater( size_t( 1));
        FragmentEnsemble split_ensemble( isolater( split_mol));
        // BCL_Debug( split_ensemble.GetSize());

        // save the final molecule(s)
        for
        (
            auto mol_itr( split_ensemble.Begin()), mol_itr_end( split_ensemble.End());
            mol_itr != mol_itr_end;
            ++mol_itr
        )
        {
          // BCL_Debug( mol_itr->GetSize());
          common_structs.PushBack( *mol_itr);
        }
      }
      return common_structs;
    }

    //! @brief reads in molecules from a given file if it is necessary
    void FragmentSplitCommonSubstructure::ReadFile() const
    {
      if( !m_FileWasRead)
      {
        io::IFStream input;
        io::File::MustOpenIFStream( input, m_File);
        m_Molecules.ReadMoreFromMdl( input);
        io::File::CloseClearFStream( input);
        m_MoleculeGraphs.AllocateMemory( m_Molecules.GetSize());
        for
        (
          FragmentEnsemble::const_iterator itr_mol( m_Molecules.Begin()), itr_mol_end( m_Molecules.End());
          itr_mol != itr_mol_end;
          ++itr_mol
        )
        {
          m_MoleculeGraphs.PushBack( m_Converter( *itr_mol));
        }
        m_FileWasRead = true;
      }
    }

    //! @brief Set the members of this property from the given LABEL
    //! @param LABEL the label to parse
    //! @param ERR_STREAM stream to write out errors to
    bool FragmentSplitCommonSubstructure::ReadInitializerSuccessHook
    (
      const util::ObjectDataLabel &LABEL,
      std::ostream &ERR_STREAM
    )
    {
      m_Converter = ConformationGraphConverter( m_AtomComparison, m_BondComparison);
      return true;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer FragmentSplitCommonSubstructure::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "splits molecules into their common substructures relative to an input set"
      );

      parameters.AddInitializer
      (
        "file",
        "file containing molecules to compare for largest common substructures",
        io::Serialization::GetAgentInputFilename( &m_File)
      );
      parameters.AddInitializer
      (
        "atom_comparison",
        "atom data that is compared to determine whether atoms are equivalent",
        io::Serialization::GetAgent( &m_AtomComparison),
        "ElementType"
      );
      parameters.AddInitializer
      (
        "bond_comparison",
        "bond data that is compared",
        io::Serialization::GetAgent( &m_BondComparison),
        "BondOrderAmideOrAromaticWithRingness"
      );
      parameters.AddInitializer
      (
        "complement",
        "return the complement of the largest common substructure instead",
        io::Serialization::GetAgent( &m_Complement),
        "0"
      );

      return parameters;
    }

  } // namespace chemistry
} // namespace bcl
