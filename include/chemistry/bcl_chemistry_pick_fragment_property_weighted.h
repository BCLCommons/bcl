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

#ifndef BCL_CHEMISTRY_PICK_FRAGMENT_PROPERTY_WEIGHTED_H_
#define BCL_CHEMISTRY_PICK_FRAGMENT_PROPERTY_WEIGHTED_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "find/bcl_find_pick_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class PickFragmentPropertyWeighted
    //! @brief Class is used for picking a fragment weighted by a stored property
    //! If the property isn't present, or isn't convertable to a float, a weight of 1.0 is used
    //! Properties are take as the absolute value (so a fragment of property -4.5 is three times as likely as 1.5)
    //!
    //! @author morettr
    //! @date 03/27/2015
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API PickFragmentPropertyWeighted :
      public find::PickInterface< const FragmentComplete &, FragmentEnsemble>
    {

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

      //! the name of the property over which the selection will be weighted.
      std::string                                            m_PropertyName;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! default constructor
      PickFragmentPropertyWeighted();

      //! @brief constructor from property name
      //! @param PROPERTY the name of the property over which to weight
      PickFragmentPropertyWeighted( const std::string &PROPERTY_NAME);

      //! virtual copy constructor
      PickFragmentPropertyWeighted *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief returns name of the property over which we're weighting
      //! @return the property name over which we're weighting
      const std::string &GetPropertyName() const 
      {
        return m_PropertyName;
      }

      //! @brief sets name of the property over which we're weighting
      //! @param PROPERTY_NAME the property name over which we're weighting
      void SetPropertyName( const std::string &PROPERTY_NAME)
      {
        m_PropertyName = PROPERTY_NAME;
      }

    ////////////////
    // operations //
    ////////////////

      //! Picks a random fragment, weighted by the appropriate property
      //! @details Fragments without the property (or without the property convertable to a float) get weight 1.0
      //! @param FRAGMENT_LIST is the list of fragments
      //! @return returns a random fragment
      const FragmentComplete &Pick( const FragmentEnsemble &FRAGMENTS) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief read from std::ostream
      //! @param OSTREAM input stream
      //! @param INDENT indentation
      //! @return ostream which was read from
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; // class PickFragmentPropertyWeighted

  } // namespace chemistry
} // namespace bcl

#endif // BCL_CHEMISTRY_PICK_FRAGMENT_PROPERTY_WEIGHTED_H_
