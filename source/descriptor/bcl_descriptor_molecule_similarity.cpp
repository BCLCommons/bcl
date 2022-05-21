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
#include "descriptor/bcl_descriptor_molecule_similarity.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_conformation_graph_converter.h"
#include "graph/bcl_graph_edge_cover_ring_perception.h"
#include "io/bcl_io_file.h"
#include "io/bcl_io_serialization.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {

    // add each of the possible instances to the enumerated instances
    const util::SiPtr< const util::ObjectInterface> MoleculeSimilarity::s_Instance
    (
      util::Enumerated< Base< chemistry::AtomConformationalInterface, float> >::AddInstance
      (
        new MoleculeSimilarity()
      )
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    // rest all histograms!
    MoleculeSimilarity::MoleculeSimilarity() :
      m_Comparer(),
      m_Molecules(),
      m_Filename()
    {
    }

    //! @brief constructor with comparer and comparison molecules
    MoleculeSimilarity::MoleculeSimilarity
    (
      const std::string &COMPARISON_TYPE,
      const chemistry::FragmentEnsemble &MOLECULES
    ) :
      m_Comparer( util::Implementation< chemistry::ConformationComparisonInterface>( COMPARISON_TYPE)),
      m_Molecules( MOLECULES),
      m_Filename( std::string())
    {
    }

    //! @brief virtual copy constructor
    MoleculeSimilarity *MoleculeSimilarity::Clone() const
    {
      return new MoleculeSimilarity( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief get name of the current class
    //! @return name of the class
    const std::string &MoleculeSimilarity::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the name of the property without any parameters
    //! @return name of the property as string
    const std::string &MoleculeSimilarity::GetAlias() const
    {
      static const std::string s_name( "MoleculeSimilarity");
      return s_name;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief calculate the descriptors
    //! @param ELEMENT the element of the sequence of interest
    //! @param STORAGE storage for the descriptor
    void MoleculeSimilarity::Calculate
    (
      linal::VectorReference< float> &STORAGE
    )
    {
      util::SiPtr< const chemistry::ConformationInterface> si_molecule( this->GetCurrentObject());
      const chemistry::ConformationInterface &molecule( *si_molecule);    //allows data type (use) of molecule

      size_t bit_position( 0);
      for
      (
        chemistry::FragmentEnsemble::const_iterator itr_mol( m_Molecules.Begin()), itr_mol_end( m_Molecules.End());
        itr_mol != itr_mol_end;
        ++itr_mol, ++bit_position
      )
      {
        STORAGE( bit_position) = ( *m_Comparer)( molecule, *itr_mol);
      }
    } // calculate

    /////////////////////////////
    /// Helper Functions/////////
    /////////////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer MoleculeSimilarity::GetSerializer() const
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
      parameters.AddInitializer
      (
        "method",
        "the method of comparison to use",
        io::Serialization::GetAgent( &m_Comparer),
        "LargestCommonSubstructureTanimoto"
      );
      return parameters;
    }

    //! @brief Set the members of this property from the given LABEL
    //! @param LABEL the label to parse
    //! @param ERR_STREAM stream to write out errors to
    bool MoleculeSimilarity::ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERR_STREAM)
    {
      io::IFStream input;
      io::File::MustOpenIFStream( input, m_Filename);
      m_Molecules = chemistry::FragmentEnsemble( input);
      io::File::CloseClearFStream( input);

      m_Comparer->PrepareEnsemble( m_Molecules);

      return true;
    }

    //! @brief hook that derived classes can override to add behavior after every time SetObject is called
    void MoleculeSimilarity::SetObjectHook()
    {
      util::SiPtr< const chemistry::ConformationInterface> si_molecule( this->GetCurrentObject());
      const chemistry::ConformationInterface &molecule( *si_molecule);
      m_Comparer->Prepare( molecule);
    }

  } // namespace descriptor
} // namespace bcl
