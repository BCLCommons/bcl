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

#ifndef BCL_DESCRIPTOR_MOLECULE_INTER_HBOND_CODE_H_
#define BCL_DESCRIPTOR_MOLECULE_INTER_HBOND_CODE_H_

// include the namespace header
#include "bcl_descriptor.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_descriptor_base_sequence.h"
#include "bcl_descriptor_molecule_3da_smooth_sign_code.h"
#include "bcl_descriptor_molecule_3da_smooth_sign_occlusion_code.h"
#include "bcl_descriptor_window.h"
#include "chemistry/bcl_chemistry_atom_conformational_interface.h"
#include "chemistry/bcl_chemistry_fragment_complete.h"
#include "util/bcl_util_implementation.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MoleculeInterHBondCode
    //! @brief code object for 3DARealSpace
    //! @details Calculates a 3DACode-style sum of products for hydrogen bond donor/acceptor pairs discretized by both
    //! distance and angle (HBD-H and HBA-H projection angle in degrees) between two interacting bodies
    //! (i.e. protein and ligand)
    //!
    //! @see @link example_descriptor_molecule_3DA_pair_convolution.cpp @endlink
    //! @author brownbp1
    //! @date Aug 04, 2020
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MoleculeInterHBondCode :
      public BaseSequence< chemistry::AtomConformationalInterface, float>
    {

    //////////
    // data //
    //////////

      // discretization data
      CheminfoProperty     m_AtomHBD;      //!< the atom property encode in the 3D autocorrelation function
      CheminfoProperty     m_AtomHBA;      //!< the atom property encode in the 3D autocorrelation function
      size_t               m_NumberSteps;  //!< number of steps in 3D autocorrelation function
      float                m_StepSize;     //!< step size for 3d autocorrelation function
      size_t                m_AngleSize;    //!< angle discretization size for 3d autocorrelation function

      //! Temporary vector used in calculations; stored here to avoid reallocating it during every call to Calculate
      linal::Vector< float> m_DiscreteCode;

      //! string to access MiscProperty in SmallMolecule
      std::string m_MiscPropertyString;

      //! Reference molecule conformer filename
      std::string m_MolBFilename;

      //! Return a file of atom indices contributing to 3da in molecule A
      std::string m_GetMolAAtomIndices;

      //! Return a file of atom indices contributing to 3da in molecule B
      std::string m_GetMolBAtomIndices;

      //! Reference molecule conformer
      chemistry::FragmentComplete m_MolB;

      linal::Vector< float> m_PocketHBD;
      linal::Vector< float> m_PocketHBA;

      //! number of features returned ( = # of features per position in the window)
      size_t m_InternalDescriptorSize;

      //! cache each protein we read in for descriptor generation
      static storage::Map
      <
        std::string,
        storage::Pair< chemistry::FragmentComplete, storage::Map< util::ObjectDataLabel, linal::Vector< float> > >
      > s_Pockets;

      //! static mutex
      static sched::Mutex s_Mutex;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      MoleculeInterHBondCode();

      //! @brief constructor from number of steps, and mapped atom property
      MoleculeInterHBondCode
      (
        const CheminfoProperty &ATOM_HBD,
        const CheminfoProperty &ATOM_HBA,
        const size_t NUMBER_STEPS = 8,
        const float STEP_SIZE = 0.50
      );

      //! @brief virtual copy constructor
      MoleculeInterHBondCode *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief get name of the current class
      //! @return name of the class
      const std::string &GetClassIdentifier() const;

      //! @brief return the data label
      //! @return data label as string
      const std::string &GetAlias() const;

      //! @brief get the feature siz  e under the normal dimension setting (e.g. GetNormalDimension())
      //! @return the feature size, assuming this feature has its normal dimension setting
      size_t GetNormalSizeOfFeatures() const
      {
        size_t angle_steps( float( 360.0) / m_AngleSize);
        return m_NumberSteps * angle_steps;
      }

      //! @brief get number of steps of code
      //! @return number of steps in 2da code
      size_t GetNumberSteps() const;

      //! @brief get step size of code
      //! @return step size of 3DA code
      float GetStepSize() const;

      //! @brief get angle size of code
      //! @return angle size of 3DA code
      size_t GetAngleSize() const;

      //! @brief get temperature of code
      //! @return const float  temperature of 3DA code
      const float &GetTemperature() const;

      //! @brief get atom property of code
      //! @return atom property mapped in 2da code
      const CheminfoProperty &GetAtomProperty() const;

    ////////////////
    // operations //
    ////////////////

    protected:

      //! @brief function to return derived-class-held implementations to this interface
      //! This allows this base class to handle mundane details like calling SetDimension and SetObject on all internal
      //! implementations
      iterate::Generic< Base< chemistry::AtomConformationalInterface, float> > GetInternalDescriptors();

      //! @brief get the cache preference under the normal dimension setting (e.g. GetType().GetDimension())
      //! @return the cache preference, assuming this feature has its normal dimension setting
      CachePreference GetNormalCachePreference() const;

      //! @brief calculate the descriptors
      //! @param STORAGE storage for the descriptor
      void Calculate( linal::VectorReference< float> &STORAGE);

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

      //! @brief Set the members of this property from the given LABEL
      //! @param LABEL the label to parse
      //! @param ERR_STREAM stream to write out errors to
      bool ReadInitializerSuccessHook
      (
        const util::ObjectDataLabel &LABEL,
        std::ostream &ERR_STREAM
      );

      //! @brief hook that derived classes can override to add behavior after every time SetObject is called
      virtual void SetObjectHook();

    }; // class MoleculeInterHBondCode

  } // namespace descriptor
} // namespace bcl

#endif // BCL_DESCRIPTOR_MOLECULE_INTER_HBOND_CODE_H_
