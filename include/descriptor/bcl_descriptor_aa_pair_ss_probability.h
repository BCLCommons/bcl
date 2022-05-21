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

#ifndef BCL_DESCRIPTOR_AA_PAIR_SS_PROBABILITY_H_
#define BCL_DESCRIPTOR_AA_PAIR_SS_PROBABILITY_H_

// include the namespace header
#include "bcl_descriptor.h"

// other forward includes
#include "model/bcl_model.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_descriptor_base_element.h"
#include "bcl_descriptor_type.h"
#include "bcl_descriptor_window_alignment_type.h"
#include "linal/bcl_linal_matrix.h"
#include "sched/bcl_sched_mutex.h"
#include "storage/bcl_storage_map.h"
#include "storage/bcl_storage_symmetric_matrix.h"
#include "storage/bcl_storage_vector.h"
#include "util/bcl_util_implementation.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AAPairSSProbability
    //! @brief Obtains probabilities of pairs of AAs from a file
    //!
    //! @see @link example_descriptor_aa_pair_ss_probability.cpp @endlink
    //! @author mendenjl
    //! @date Aug 27, 2013
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API AAPairSSProbability :
      public BaseElement< biol::AABase, float>
    {
    public:

    //////////
    // data //
    //////////

      //! custom SSE info type
      enum SSInfoType
      {
        e_Helix,         //!< P(central residue is in a helix)
        e_Strand,        //!< P(central residue is in a strand)
        e_Coil,          //!< P(central residue is in a coil)
        e_ToHelix,       //!< P(SS transition to helix over the stencil)
        e_ToStrand,      //!< P(SS transition to strand over the stencil)
        e_ToCoil,        //!< P(SS transition to coil over the stencil)
        e_FromHelix,     //!< P(SS transition from helix over the stencil)
        e_FromStrand,    //!< P(SS transition from strand over the stencil)
        e_FromCoil,      //!< P(SS transition from coil over the stencil)
        s_NumberSSInfoTypes
      };

      //! @brief SSInfoType as string
      //! @param SS_INFO_TYPE the type of interest
      //! @return the SSInfoType as string
      static const std::string &GetSSInfoTypeString( const SSInfoType &SS_INFO_TYPE);

      //! @brief SSInfoType directionality (relative to the central residue that the descriptor considers)
      //! @param SS_INFO_TYPE the type of interest
      //! @return -1 if the left is considered, 1 if the right is considered, 0 if both are considered
      static int GetDirectionality( const SSInfoType &SS_INFO_TYPE);

      //! @brief Get the column offset in the input file for the given type
      //! @param SS_INFO_TYPE the type of interest
      //! @return Column offset (4 for helical types, 7 for strand types, 10 for coil types)
      static size_t GetColumnOffset( const SSInfoType &SS_INFO_TYPE);

      //! @brief get a vector that takes as input any letter, and returns the corresponding AA type index
      //! @return the vector
      static const storage::Vector< size_t> &GetAATypeId();

      //! SSInfoTypeEnum simplifies the usage of the SSInfoType enum of this class
      typedef util::WrapperEnum< SSInfoType, &GetSSInfoTypeString, s_NumberSSInfoTypes> SSInfoTypeEnum;

    private:

    //////////
    // data //
    //////////

      //! Stencil of residues (does not include the central residue)
      //! e.g. 1,4 means that probabilities will be calculated
      //! P(In Helix | Res at -4 Res at 0) * P(In Helix | Res at -1 Res at 0) * P(In Helix | Res at 4 Res at 0) * P(In Helix | Res at 1 Res at 0)
      linal::Vector< size_t> m_Stencil;

      //! Actual type of information calculated
      SSInfoTypeEnum m_Type;

      //! Actual file name for the statistics file
      std::string m_Filename;

      WindowAlignmentEnum m_ForceAlignment; //!< Force the type to use left or right alignment even if m_Type says otherwise

      //! descriptor used to calculate AA type at each point in the stencil; must return 20 values
      util::Implementation< Base< biol::AABase, float> > m_TypeCalculator;

      //! Space to store output of m_TypeCalculator on each stenciled location
      storage::Vector< linal::Vector< float> >    m_TypeCalculatorOutput;
      storage::Vector< storage::Vector< size_t> > m_SignificantTypes;

      //! Matrix of probabilities at each distance for the given sstype
      util::SiPtr< const storage::Vector< linal::Matrix< float> > > m_Probabilities;

      //! Mutex for access to s_ProbabilityMapStorage
      static sched::Mutex s_ProbabilityMapMutex;

      //! static map holding all the probabilities from each file; key is label for this object
      static storage::Map< util::ObjectDataLabel, storage::Vector< linal::Matrix< float> > > s_ProbabilityMapStorage;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      AAPairSSProbability() :
        m_ForceAlignment( s_NumberOfWindowAlignments)
      {
      }

      //! @brief virtual copy constructor
      AAPairSSProbability *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief return the data label
      //! @return data label as string
      const std::string &GetAlias() const;

    //////////////////////
    // helper functions //
    //////////////////////

    protected:

      //! @brief get the feature size under the normal dimension setting (e.g. GetNormalDimension())
      //! @return the feature size, assuming this feature has its normal dimension setting
      size_t GetNormalSizeOfFeatures() const;

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    private:

      //! @brief calculate the descriptors
      //! @param ELEMENT: the element pair of interest
      //! @param STORAGE storage for the descriptor
      //! @return true, if the calculation was successful
      void Calculate
      (
        const iterate::Generic< const biol::AABase> &ELEMENT,
        linal::VectorReference< float> &STORAGE
      );

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief Set the members of this property from the given LABEL
      //! @param LABEL the label to parse
      //! @param ERR_STREAM stream to write out errors to
      bool ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERR_STREAM);

      //! @brief function to return derived-class-held implementations to this interface
      //! This allows this base class to handle mundane details like calling SetDimension and SetObject on all internal
      //! implementations
      iterate::Generic< Base< biol::AABase, float> > GetInternalDescriptors()
      {
        return iterate::Generic< Base< biol::AABase, float> >( &m_TypeCalculator, &m_TypeCalculator + 1);
      }

      //! @brief get the cache preference under the normal dimension setting (e.g. GetType().GetDimension())
      //! @return the cache preference, assuming this feature has its normal dimension setting
      CachePreference GetNormalCachePreference() const
      {
        return e_PreferCache;
      }

      //! @brief get the effective directionality (considering implicit direction of m_Type, and m_ForceAlignment)
      //! @return -1 for left, 0 for center, 1 for right
      int GetDirectionality() const;

    }; // class AAPairSSProbability

  } // namespace descriptor
} // namespace bcl

#endif // BCL_DESCRIPTOR_AA_PAIR_SS_PROBABILITY_H_
