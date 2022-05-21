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
#include "descriptor/bcl_descriptor_reaction_structure_search.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_atom_constitutional_interface.h"
#include "chemistry/bcl_chemistry_fragment_ensemble.h"
#include "graph/bcl_graph_edge_cover_ring_perception.h"
#include "io/bcl_io_file.h"
#include "io/bcl_io_serialization.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {

    // add each of the possible instances to the enumerated instances
    const util::SiPtr< const util::ObjectInterface> ReactionStructureSearch::s_Instance
    (
      util::Enumerated< Base< chemistry::AtomConformationalInterface, float> >::AddInstance
      (
        new ReactionStructureSearch()
      )
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    // rest all histograms!
    ReactionStructureSearch::ReactionStructureSearch()
    {
    }

    //! @brief virtual copy constructor
    ReactionStructureSearch *ReactionStructureSearch::Clone() const
    {
      return new ReactionStructureSearch( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief get name of the current class
    //! @return name of the class
    const std::string &ReactionStructureSearch::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the name of the property without any parameters
    //! @return name of the property as string
    const std::string &ReactionStructureSearch::GetAlias() const
    {
      static const std::string s_name( "ReactionStructureSearch");
      return s_name;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief calculate the descriptors
    //! @param ELEMENT the element of the sequence of interest
    //! @param STORAGE storage for the descriptor
    void ReactionStructureSearch::Calculate
    (
      linal::VectorReference< float> &STORAGE
    )
    {
      util::SiPtr< const chemistry::ConformationInterface> si_molecule( this->GetCurrentObject());
      const chemistry::ConformationInterface &molecule( *si_molecule);    //allows data type (use) of molecule

      chemistry::ReactionStructure mol_struct( molecule);

      size_t bit_position( 0);
      for
      (
        storage::Vector< chemistry::ReactionStructure>::iterator
        itr_frags( m_Structs.Begin()), itr_frags_end( m_Structs.End());
        itr_frags != itr_frags_end;
        ++itr_frags, ++bit_position
      )
      {
        //check for isomorphism
        if( mol_struct.Contains( *itr_frags))
        {
          STORAGE( bit_position) = 1;
        }
      }
    } // calculate

    /////////////////////////////
    /// Helper Functions/////////
    /////////////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer ReactionStructureSearch::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Searches for specified substructures within the query molecule. Returns 1 if a structure exists in the query "
        "molecule"
      );
      parameters.AddInitializer
      (
        "filename",
        "File containing molecules to compare with each query molecule",
        io::Serialization::GetAgentInputFilename( &m_Filename)
      );
      return parameters;
    }

    //! @brief Set the members of this property from the given LABEL
    //! @param LABEL the label to parse
    //! @param ERR_STREAM stream to write out errors to
    bool ReactionStructureSearch::ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERR_STREAM)
    {
      io::IFStream input;
      io::File::MustOpenIFStream( input, m_Filename);
      chemistry::FragmentEnsemble fragments( input);
      io::File::CloseClearFStream( input);

      m_Structs.AllocateMemory( fragments.GetSize());
      for
      (
        chemistry::FragmentEnsemble::const_iterator itr( fragments.Begin()), itr_end( fragments.End());
        itr != itr_end;
        ++itr
      )
      {
        m_Structs.PushBack( chemistry::ReactionStructure( *itr));
      }
      return true;
    }

  } // namespace descriptor
} // namespace bcl
