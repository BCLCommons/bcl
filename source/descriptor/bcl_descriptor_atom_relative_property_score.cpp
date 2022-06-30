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
#include "chemistry/bcl_chemistry_configurational_bond_type_data.h"
#include "chemistry/bcl_chemistry_molecule_feature_mapper.h"
#include "io/bcl_io_file.h"
#include "io/bcl_io_ifstream.h"
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "descriptor/bcl_descriptor_atom_relative_property_score.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_conformation_graph_converter.h"
#include "graph/bcl_graph_connectivity.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    AtomRelativePropertyScore::AtomRelativePropertyScore() :
        m_Descriptor( "XLogP"),
        m_ReferenceMolsFilename( ""),
        m_ReferenceMols( chemistry::FragmentEnsemble())
    {
    }

    //! @brief construct with property string
    AtomRelativePropertyScore::AtomRelativePropertyScore( const CheminfoProperty &PROPERTY) :
        m_Descriptor( PROPERTY),
        m_ReferenceMolsFilename( ""),
        m_ReferenceMols( chemistry::FragmentEnsemble())
    {
    }

    //! @brief construct with reference filename
    AtomRelativePropertyScore::AtomRelativePropertyScore( const std::string &REFERENCE_FILENAME) :
        m_Descriptor( "XLogP"),
        m_ReferenceMolsFilename( REFERENCE_FILENAME),
        m_ReferenceMols( chemistry::FragmentEnsemble())
    {
    }

    //! @brief virtual copy constructor
    //! @return pointer to new AtomTopologicalPolarSurface
    AtomRelativePropertyScore *AtomRelativePropertyScore::Clone() const
    {
      return new AtomRelativePropertyScore( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief get name of the current class
    //! @return name of the class
    const std::string &AtomRelativePropertyScore::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the name of the property without any parameters
    //! @return name of the property as string
    const std::string &AtomRelativePropertyScore::GetAlias() const
    {
      static const std::string s_name( "Atom_RelativePropertyScore");
      return s_name;
    }

  ////////////////
  // operations //
  ////////////////

    void AtomRelativePropertyScore::SetReferenceMols( chemistry::FragmentEnsemble REFERNCE_MOLS)
    {
      m_ReferenceMols = REFERNCE_MOLS;
    }

    //! @brief calculate the descriptors
    //! @param ELEMENT the element of the sequence of interest
    //! @param STORAGE storage for the descriptor
    void AtomRelativePropertyScore::Calculate
    (
      const iterate::Generic< const chemistry::AtomConformationalInterface> &ELEMENT,
      linal::VectorReference< float> &STORAGE
    )
    {
      // Get our molecule of interest
      util::SiPtr< const chemistry::ConformationInterface> conformation( GetCurrentObject());

      // build feature mapper object with the scoring object (e.g. AffinityNet)
      std::map< size_t, float>::iterator features_itr( m_FeatureMap[ 0].find( conformation->GetAtomIndex( *ELEMENT)));

      // set atom dddProperty value
      STORAGE( 0) = features_itr->second;
    } // Recalculate

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief Set the members of this property from the given LABEL
    //! @param LABEL the label to parse
    //! @param ERR_STREAM stream to write out errors to
    bool AtomRelativePropertyScore::ReadInitializerSuccessHook
    (
      const util::ObjectDataLabel &LABEL,
      std::ostream &ERR_STREAM
    )
    {
      return true;
    }

    //! @brief hook that derived classes can override to add behavior after every time SetObject is called
    void AtomRelativePropertyScore::SetObjectHook()
    {
      //      m_Descriptor = descriptor::CheminfoProperty();
      if( m_ReferenceMolsFilename.size())
      {
        io::IFStream input;
        io::File::MustOpenIFStream( input, m_ReferenceMolsFilename);
        m_ReferenceMols.ReadMoreFromMdl( input);
      }

      // generate feature map for molecule
      util::SiPtr< const chemistry::ConformationInterface> conformation( GetCurrentObject());

      static chemistry::MoleculeFeatureMapper feature_mapper( m_Descriptor.GetLabel());
      m_FeatureMap =
          feature_mapper.CompareSubstructuresRigorous
          (
            *conformation,
            m_ReferenceMols,
            chemistry::ConformationGraphConverter::e_ElementType,
            chemistry::ConfigurationalBondTypeData::e_BondOrderAmideOrAromaticWithRingness,
            false,
            true,
            false,
            util::GetUndefinedSize_t(),
            false
          );
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer AtomRelativePropertyScore::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription( "Compute per-atom contribution of QSAR score relative to a provided scaffold");
      parameters.AddInitializer
      (
        "property",
        "the property to use for scoring",
        io::Serialization::GetAgent( &m_Descriptor),
        "XLogP"
      );
      parameters.AddInitializer
      (
        "reference_mols",
        "the reference molecule(s) against which target molecules will be scored",
        io::Serialization::GetAgent( &m_ReferenceMolsFilename),
        ""
      );
      return parameters;
    }

  } // namespace descriptor
} // namespace bcl
