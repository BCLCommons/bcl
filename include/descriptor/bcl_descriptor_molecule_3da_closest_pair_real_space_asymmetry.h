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

#ifndef BCL_DESCRIPTOR_MOLECULE_3DA_CLOSEST_PAIR_REAL_SPACE_ASYMMETRY_H_
#define BCL_DESCRIPTOR_MOLECULE_3DA_CLOSEST_PAIR_REAL_SPACE_ASYMMETRY_H_

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
    //! @class Molecule3DAClosestPairRealSpaceAsymmetry
    //! @brief code object for 3DARealSpace
    //! @details Calculates the 3DASign of the contacting residues between a protein and a ligand
    //!  such that the 3DA distance bins are partitioned across contact distance.
    //!
    //! @see @link example_descriptor_molecule_3DA_pair_convolution.cpp @endlink
    //! @author brownbp1, mendenjl
    //! @date Sep 21, 2019
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API Molecule3DAClosestPairRealSpaceAsymmetry :
      public BaseSequence< chemistry::AtomConformationalInterface, float>
    {

    //////////
    // data //
    //////////

      CheminfoProperty     m_AtomPropertyA; //!< the atom property encoded in the 3D autocorrelation function of the small molecule
      CheminfoProperty     m_AtomPropertyB; //!< the atom property encoded in the 3D autocorrelation function of the protein pocket
      size_t               m_NumberSteps;  //!< number of steps in 3D autocorrelation function
      float                m_StepSize;     //!< step size for 3d autocorrelation function
      float                m_Temperature;  //!< Exponent of the gaussian function
      bool                 m_Smooth;       //!< Whether to perform smoothing; otherwise, linear kernel is used
      bool                 m_Interpolate;  //!< Whether to perform interpolation; otherwise, the nearest bin will always be selected

      //! Temporary vector used in calculations; stored here to avoid reallocating it during every call to Calculate
      linal::Vector< float> m_DiscreteCode;

      //! Cached vector of smoothing coefficients; used to avoid excessive computations of the gaussian kernel
      linal::VectorConstReference< float> m_SmoothingCoefficients;

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

      linal::Vector< float> m_PocketProperties;

      //! number of features returned ( = # of features per position in the window)
      size_t m_InternalDescriptorSize;

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
      Molecule3DAClosestPairRealSpaceAsymmetry();

      //! @brief constructor from number of steps, and mapped atom property
      Molecule3DAClosestPairRealSpaceAsymmetry
      (
        const CheminfoProperty &ATOM_PROPERTY_A,
        const CheminfoProperty &ATOM_PROPERTY_B,
        const size_t NUMBER_STEPS = 14,
        const float STEP_SIZE = 0.50,
        const float TEMPERATURE = 5.0,
        const bool SMOOTH = true,
        const bool INTERPOLATE = true
      );

      //! @brief virtual copy constructor
      Molecule3DAClosestPairRealSpaceAsymmetry *Clone() const;

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
        return m_NumberSteps * 4;
      }

      //! @brief get number of steps of code
      //! @return number of steps in 2da code
      size_t GetNumberSteps() const;

      //! @brief get step size of code
      //! @return step size of 3DA code
      float GetStepSize() const;

      //! @brief get temperature of code
      //! @return const float  temperature of 3DA code
      const float &GetTemperature() const;

      //! @brief get atom property A of code
      //! @return atom property A mapped in 3da code
      const CheminfoProperty &GetAtomPropertyA() const;

      //! @brief get atom property A of code
      //! @return atom property A mapped in 3da code
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

      //! @brief create a gaussian-smoothed signal from m_DiscreteCode and store it in STORAGE
      //! @param STORAGE storage for the gaussian-smoothed signal
      void Smooth( linal::VectorReference< float> &STORAGE) const;

      //! @brief add an observed distance/property value to m_DiscreteCode
      //! @param DISTANCE actual distance of the two atoms
      //! @param PROP_A property from atom A
      //! @param PROP_B property from atom B
      void Accumulate
      (
        const float &DISTANCE,
        const float &PROP_A,
        const float &PROP_B
      );

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

    }; // class Molecule3DAClosestPairRealSpaceAsymmetry

  } // namespace descriptor
} // namespace bcl

#endif // BCL_DESCRIPTOR_MOLECULE_3DA_CLOSEST_PAIR_REAL_SPACE_ASYMMETRY_H_
