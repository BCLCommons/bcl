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
//BCL_StaticInitializationFiascoFinder

// include header of this class
#include "descriptor/bcl_descriptor_umol2d.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_atom_environment_bender.h"
#include "chemistry/bcl_chemistry_molecule_environment.h"
#include "io/bcl_io_file.h"
#include "io/bcl_io_serialization.h"
#include "score/bcl_score.h"
#include "util/bcl_util_string_functions.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {

  //////////
  // data //
  //////////

    // add the interface to the set of known implementations
    const util::SiPtr< const util::ObjectInterface> UMol2D::s_Instance
    (
      util::Enumerated< Base< chemistry::AtomConformationalInterface, float> >::AddInstance( new UMol2D())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    UMol2D::UMol2D() : m_AtomType( chemistry::AtomEnvironmentBender::e_Atom), m_AeNumber( 574), m_Sphere( 1)
    {
    }

    //! @brief virtual copy constructor
    //! @return pointer to new MoleculeRotatableBonds
    UMol2D *UMol2D::Clone() const
    {
      return new UMol2D( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief get name of the current class
    //! @return name of the class
    const std::string &UMol2D::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the name of the property without any parameters
    //! @return name of the property as string
    const std::string &UMol2D::GetAlias() const
    {
      static const std::string s_name( "UMol2D");
      return s_name;
    }

    //! @brief return the name of the property without any parameters
    //! @return name of the property as string
    storage::Vector< chemistry::AtomEnvironmentBender> UMol2D::GetAEs( const t_AtomTypeEnum &ATOM_TYPE, const size_t &SPHERE)
    {
      //! read in the list of AE string forms to convert to AEs
      io::IFStream read;
      io::File::MustOpenIFStream
      (
        read, score::Score::AddHistogramPath
        (
          ATOM_TYPE.GetString() + "_AE_" + util::Format()( SPHERE) + "_bond_sorted.txt"
        )
      );

      //! read lines by lines and store all the hashed strings of AEs into AE_strings
      storage::Vector< std::string> ae_strings( util::StringLineListFromIStream( read));
      io::File::CloseClearFStream( read);
      storage::Vector< chemistry::AtomEnvironmentBender> AEs;

      // only allocate m_aeNumber to AEs when it is smaller than 1000
      AEs.AllocateMemory( ae_strings.GetSize());

      //! generate AEs from their string representations and add into AEs vector
      for
      (
        auto string_iter( ae_strings.Begin()), string_iter_end( ae_strings.End());
        string_iter != string_iter_end; ++string_iter
      )
      {
        AEs.PushBack( chemistry::AtomEnvironmentBender( ATOM_TYPE, *string_iter));
      }
      return AEs;
    }

    //! @brief get the feature size under the normal dimension setting (e.g. GetNormalDimension())
    //! @return the feature size, assuming this feature has its normal dimension setting
    size_t UMol2D::GetNormalSizeOfFeatures() const
    {
      return m_AeNumber;
    }

    //! @brief get the cache preference under the normal dimension setting (e.g. GetType().GetDimension())
    //! @return the cache preference, assuming this feature has its normal dimension setting
    CachePreference UMol2D::GetNormalCachePreference() const
    {
      return e_PreferCache;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief calculate the descriptors
    //! @param STORAGE storage for the descriptor
    //! @return update the STORAGE to record the counts of each AE in the AE set of the molecule
    void UMol2D::Calculate( linal::VectorReference< float> &STORAGE)
    {
      std::fill( STORAGE.Begin(), STORAGE.End(), 0.0);
      util::SiPtr< const chemistry::ConformationInterface> si_molecule( this->GetCurrentObject());
      const chemistry::MoleculeEnvironment molecule( m_AtomType, *si_molecule);
      //s_AE_common_set is a vector of the selected most common atom environment of 2 different atom types
      static storage::Vector< chemistry::AtomEnvironmentBender> s_AE_common_set[ 2] =
      {
        GetAEs( chemistry::AtomEnvironmentBender::e_Element, m_Sphere),
        GetAEs( chemistry::AtomEnvironmentBender::e_Atom, m_Sphere)
      };
      //!> loop through the molecular environment
      for
      (
        storage::Vector< chemistry::AtomEnvironmentBender>::const_iterator itr_env( molecule.GetMoleculeEnvironment().Begin());
        itr_env != molecule.GetMoleculeEnvironment().End(); ++itr_env
      )
      {
        //!> compare each AE of the molecule to all AEs in the AEs set
        linal::VectorReference< float>::iterator out_iter( STORAGE.Begin());
        size_t count( 0);
        for
        (
          auto itr( s_AE_common_set[ m_AtomType].Begin()); count < m_AeNumber;
          ++count, ++itr, ++out_iter
        )
        {
          // update the count of each atom environment accordingly
          if( ( *itr_env) == ( *itr))
          {
            ++( *out_iter);
          }
        }
      }
    } // Recalculate

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer UMol2D::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription( "calculates the number of certain types of atom environments");
      parameters.AddInitializer
      (
        "atom hashing type",
        "Choose one of two: Element, Atom. Atom by default",
        io::Serialization::GetAgent( &m_AtomType),
        "Atom"
      );
      parameters.AddInitializer
      (
        "feature size",
        "get the feature size under the normal dimension setting. "
        "Range: 1-250 for Element(height=1), 1-574 for Atom(height=1), 1-5117 for Element(height=2), 1-8080 for Atom(height=2)"
        "Default: 574",
        io::Serialization::GetAgent( &m_AeNumber),
        "574"
      );
      parameters.AddInitializer
      (
        "Atom environment height",
        "Number of bond sphere from the center atom "
        "Range: 1 or 2"
        "Default: 1",
        io::Serialization::GetAgent( &m_Sphere),
        "1"
      );
      return parameters;
    }

  } // namespace descriptor
} // namespace bcl
