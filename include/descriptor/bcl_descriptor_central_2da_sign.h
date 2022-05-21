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

#ifndef BCL_DESCRIPTOR_CENTRAL_2DA_SIGN_H_
#define BCL_DESCRIPTOR_CENTRAL_2DA_SIGN_H_

// include the namespace header
#include "bcl_descriptor.h"
#include "bcl_descriptor_base_sequence.h"

// include other forward headers - sorted alphabetically
#include "chemistry/bcl_chemistry.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_descriptor_base_element.h"
#include "chemistry/bcl_chemistry_atom_conformational_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Central2DASign
    //! @brief Measure 2DASign of property at variable distances from molecule topological center
    //!
    //! @see @link example_descriptor_central_2da_sign.cpp @endlink
    //! @author brownbp1, mendenjl
    //! @date Jun 21, 2020
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API Central2DASign :
      public BaseSequence< chemistry::AtomConformationalInterface, float>
    {

    //////////
    // data //
    //////////

    private:

      // data for the outer array
      CheminfoProperty                    m_AtomProperty;               //!< the atom property encode in the 3D autocorrelation function
      size_t                              m_MaxBondDistanceFromCenter;  //!< number of steps in

      // data for the signed 2DA
      size_t                              m_Number2DASteps;             //!< number of steps in 2D autocorrelation function
      float                               m_Temperature;                //!< Exponent of the Gaussian function
      bool                                m_Smooth;                     //!< Whether to perform smoothing; otherwise, linear kernel is used
      Molecule2DASmoothSignCode           m_2DASmoothSign;              //! Used to compute the 2DA

      graph::ConstGraph< size_t, size_t>  m_Graph;                      //!< Graph of the current molecule

      //bool                                m_Smooth;                     //!< Whether to perform smoothing; otherwise, linear kernel is used
      //Molecule2DASmoothSignCode           m_2DASmoothSign;              //! Used to compute the 2DA

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      Central2DASign();

      //! @brief constructor from number of steps, and mapped atom property
      Central2DASign
      (
        const CheminfoProperty &ATOM_PROPERTY,
        const size_t MAX_CENTER_BOND_DISTANCE,
        const size_t NUMBER_2DA_STEPS,
        const float &TEMPERATURE,
        const bool &SMOOTH
      );

      //! @brief virtual copy constructor
      Central2DASign *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief get name of the current class
      //! @return name of the class
      const std::string &GetClassIdentifier() const;

      //! @brief return the data label
      //! @return data label as string
      const std::string &GetAlias() const;

      //! @brief get the feature size under the normal dimension setting (e.g. GetNormalDimension())
      //! @return the feature size, assuming this feature has its normal dimension setting
      size_t GetNormalSizeOfFeatures() const
      {
        return m_Number2DASteps * 3 * ( m_MaxBondDistanceFromCenter + 1);
      }

    ////////////////
    // operations //
    ////////////////

    protected:

      //! @brief function to return derived-class-held implementations to this interface
      //! This allows this base class to handle mundane details like calling SetDimension and SetObject on all internal
      //! implementations
      iterate::Generic< Base< chemistry::AtomConformationalInterface, float> > GetInternalDescriptors();

      //! @brief calculate the descriptors
      //! @param ELEMENT the element of the sequence of interest
      //! @param STORAGE storage for the descriptor
      void Calculate
      (
        linal::VectorReference< float> &STORAGE
      );

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief return the bond girth from every atom
      //! @return bond girth at each atom index
      storage::Vector< size_t> CalculateBondGirths( const chemistry::FragmentComplete &MOLECULE);

      //! @brief return the minimum bond girth
      //! @return shortest bond girth
      size_t CalculateMinBondGirth( const chemistry::FragmentComplete &MOLECULE);

      //! @brief Set the members of this property from the given LABEL
      //! @param LABEL the label to parse
      //! @param ERR_STREAM stream to write out errors to
      bool ReadInitializerSuccessHook
      (
        const util::ObjectDataLabel &LABEL,
        std::ostream &ERR_STREAM
      );

      //! @brief hook that derived classes can override to add behavior after every time SetObject is called
      void SetObjectHook();

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

      //! @brief get the cache preference under the normal dimension setting (e.g. GetType().GetDimension())
      //! @return the cache preference, assuming this feature has its normal dimension setting
      CachePreference GetNormalCachePreference() const
      {
        return e_PreferCache;
      }

    }; // class Central2DASign

  } // namespace descriptor
} // namespace bcl

#endif // BCL_DESCRIPTOR_CENTRAL_2DA_SIGN_H_
