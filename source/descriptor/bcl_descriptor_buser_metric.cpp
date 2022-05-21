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

// include header of this class
#include "descriptor/bcl_descriptor_buser_metric.h"
// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_fragment_ensemble.h"
#include "descriptor/bcl_descriptor_base.h"
#include "io/bcl_io_file.h"
#include "io/bcl_io_serialization.h"
#include "score/bcl_score.h"

namespace bcl
{
  namespace descriptor
  {

  //////////
  // data //
  //////////

    // add the interface to the set of known implementations
    const util::SiPtr< const util::ObjectInterface> BuserMetric::s_Instance
    (
      util::Enumerated< Base< chemistry::AtomConformationalInterface, float> >::AddInstance( new BuserMetric())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    BuserMetric::BuserMetric() : m_AtomType( chemistry::AtomEnvironmentBender::e_Element), m_ActiveNum( 1)
    {
    }

    //! @brief virtual copy constructor
    //! @return pointer to new MoleculeRotatableBonds
    BuserMetric *BuserMetric::Clone() const
    {
      return new BuserMetric( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief get name of the current class
    //! @return name of the class
    const std::string &BuserMetric::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the name of the property without any parameters
    //! @return name of the property as string
    const std::string &BuserMetric::GetAlias() const
    {
      static const std::string s_name( "BuserMetric");
      return s_name;
    }

    //! @brief get the feature size under the normal dimension setting (e.g. GetNormalDimension())
    //! @return the feature size, assuming this feature has its normal dimension setting
    size_t BuserMetric::GetNormalSizeOfFeatures() const
    {
      return m_ActiveNum;
    }

    //! @brief get the cache preference under the normal dimension setting (e.g. GetType().GetDimension())
    //! @return the cache preference, assuming this feature has its normal dimension setting
    CachePreference BuserMetric::GetNormalCachePreference() const
    {
      return e_PreferCache;
    }
  ///////////////
  // operations //
  ////////////////

    //! @brief calculate the descriptors
    //! @param STORAGE storage for the descriptor
    //! @return update the STORAGE to record the counts of each AE in the AE set of the molecule
    void BuserMetric::Calculate( linal::VectorReference< float> &STORAGE)
    {
      util::SiPtr< const chemistry::ConformationInterface> si_molecule( this->GetCurrentObject());

      //! variable to check if the molecule is one of the actives
      //bool overlapped(false);
      size_t overlap_index( m_ActiveNum);

      // Calculates the Buser score between this molecule and the active molecule
      chemistry::MoleculeEnvironment molecule( m_AtomType, *si_molecule);
      //BCL_Debug(m_Actives);
      auto out_iter( STORAGE.Begin());
      for( auto iter( m_Actives.Begin()); iter != m_Actives.End(); ++iter, ++out_iter)
      {
        if( *iter != molecule)
        {
          *out_iter = float( molecule.BuserScore( *iter));
        }
        //! when the molecule is one of the actives
        else
        {
          *out_iter = float( 0.0);
          overlap_index = std::distance( STORAGE.Begin(), out_iter);
        }
      }

      //! if the molecule is one of the actives
      if( overlap_index < m_ActiveNum && m_ActiveNum > 1)
      {
        //! set buser value for the overlapped active to be the average of all the other actives
        float average( STORAGE.Sum() / ( m_ActiveNum - 1));
        STORAGE( overlap_index) = average;
      }

    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer BuserMetric::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription( "calculates the buser similarity score between ");
      parameters.AddInitializer
      (
        "atom hashing type",
        "Choose one of four: Element, ElemRC, Atom, AtomRC. element by default",
        io::Serialization::GetAgent( &m_AtomType),
        "Element"
      );
      parameters.AddInitializer
      (
        "filename",
        "File containing active compounds to compare with each query molecule",
        io::Serialization::GetAgentInputFilename( &m_FileName)
      );

      return parameters;
    }

    //! @brief Set the members of this property from the given LABEL
    //! @param LABEL the label to parse
    //! @param ERR_STREAM stream to write out errors to
    bool BuserMetric::ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERR_STREAM)
    {

      //! read the input file containing active molecules
      io::IFStream input;
      io::File::MustOpenIFStream( input, m_FileName);
      chemistry::FragmentEnsemble fragments( input);
      io::File::CloseClearFStream( input);
      m_Actives.AllocateMemory( fragments.GetSize());
      // convert to molecule environments and add them into m_ActiveNum
      for
      (
          chemistry::FragmentEnsemble::const_iterator itr( fragments.Begin()), itr_end( fragments.End());
          itr != itr_end;
          ++itr
      )
      {
        m_Actives.PushBack( chemistry::MoleculeEnvironment( m_AtomType, *itr));
      }
      m_ActiveNum = m_Actives.GetSize();
      return true;
    }
  } // namespace descriptor
} // namespace bcl
