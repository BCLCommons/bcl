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

#ifndef BCL_DESCRIPTOR_MOLECULE_3DA_PAIR_CONVOLUTION_ASYMMETRY_H_
#define BCL_DESCRIPTOR_MOLECULE_3DA_PAIR_CONVOLUTION_ASYMMETRY_H_

// include the namespace header
#include "bcl_descriptor.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_descriptor_base_sequence.h"
#include "bcl_descriptor_cache_preference.h"
#include "bcl_descriptor_molecule_3da_smooth_sign_code.h"
#include "bcl_descriptor_molecule_3da_smooth_sign_occlusion_code.h"
#include "bcl_descriptor_window.h"
#include "bcl_descriptor_window_weighting_interface.h"
#include "chemistry/bcl_chemistry_atom_conformational_interface.h"
#include "chemistry/bcl_chemistry_fragment_complete.h"
#include "sched/bcl_sched_mutex.h"
#include "util/bcl_util_implementation.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Molecule3DAPairConvolutionAsymmetry
    //! @brief code object for 3DAPairConvolution
    //! @details Relates the 3DASign of two molecules by correlating matching signed distance bins.
    //!          The 3DASign of the first molecule is convolved prior to multiplication with the second molecule.
    //!
    //! @see @link example_descriptor_molecule_3DA_pair_convolution.cpp @endlink
    //! @author brownbp1, mendenjl
    //! @date Oct 10, 2019
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API Molecule3DAPairConvolutionAsymmetry :
      public BaseSequence< chemistry::AtomConformationalInterface, float>
    {

    //////////
    // data //
    //////////

      //!< size (in multiples of StepSize) of the window around each 3da bin
      size_t m_WindowSize;

      //! string to access MiscProperty in SmallMolecule
      std::string m_MiscPropertyString;

      //! string naming pocket
      std::string m_PocketName;

      //! Reference molecule conformer filename
      std::string m_MolBFilename;

      //! Initial 3DA for molecule 2 property to be used in convolution
      linal::Vector< float> m_MolB3DA;

      //! Used to compute the individual 3DAs of our molecules and our reference molecule conformer
      Molecule3DASmoothSignCode m_3DASmoothSignA;
      Molecule3DASmoothSignCode m_3DASmoothSignB;

      //! object used to create the window weights
      util::Implementation< WindowWeightingInterface> m_WindowWeightsCreator;

      //! actual weights for the window
      linal::Vector< float> m_WindowWeights;

    protected:

      //! Reference molecule conformer
      chemistry::FragmentComplete m_MolB;

      //! number of features returned ( = # of features per position in the window)
      size_t m_InternalDescriptorSize;

      //! the atom property encode in the 3D autocorrelation function of the ligand
      CheminfoProperty m_AtomPropertyA;

      //! the atom property encode in the 3D autocorrelation function of the pocket
      CheminfoProperty m_AtomPropertyB;

      //!< number of steps in 3D autocorrelation function
      size_t m_NumberSteps;

      //!< step size for 3d autocorrelation function
      float m_StepSize;

      //! cache each protein we read in for descriptor generation
      static storage::Map< std::string, storage::Pair< chemistry::FragmentComplete, storage::Map< util::ObjectDataLabel, linal::Vector< float> > > > s_Pockets;

      //! static mutex
      static sched::Mutex s_Mutex;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      Molecule3DAPairConvolutionAsymmetry();

      //! @brief constructor from number of steps, and mapped atom property
      Molecule3DAPairConvolutionAsymmetry
      (
        const CheminfoProperty &ATOM_PROPERTY_A,
        const CheminfoProperty &ATOM_PROPERTY_B,
        const size_t NUMBER_STEPS = 20,
        const float STEP_SIZE = 0.50,
        const size_t WINDOW_SIZE = 4.0
      );

      //! @brief virtual copy constructor
      Molecule3DAPairConvolutionAsymmetry *Clone() const;

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
          return m_3DASmoothSignA.GetNormalSizeOfFeatures() * 3;
      }

      //! @brief get number of steps of code
      //! @return number of steps in 2da code
      size_t GetNumberSteps() const;

      //! @brief get step size of code
      //! @return step size of 3DA code
      float GetStepSize() const;

      //! @brief get atom property of code
      //! @return atom property mapped in 3da code
      const CheminfoProperty &GetAtomPropertyA() const;

      //! @brief get atom property of code
      //! @return atom property mapped in 3da code
      const CheminfoProperty &GetAtomPropertyB() const;

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

      //! @brief reduce pocket to atoms contacting the ligand
      virtual void ReduceProteinPocket( chemistry::FragmentComplete &POCKET);

    }; // class Molecule3DAPairConvolutionAsymmetry

  } // namespace descriptor
} // namespace bcl

#endif // BCL_DESCRIPTOR_MOLECULE_3DA_PAIR_CONVOLUTION_ASYMMETRY_H_
