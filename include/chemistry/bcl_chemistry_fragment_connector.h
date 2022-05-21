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

#ifndef BCL_CHEMISTRY_FRAGMENT_CONNECTOR_H_
#define BCL_CHEMISTRY_FRAGMENT_CONNECTOR_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_fragment_complete.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class FragmentConnector
    //! @brief class creates intermediate fully connected fragments during molecule assembly by joining fragments at given atoms
    //! @details given a fragment to be attached to currently held held fragment, this class attaches the incoming fragment
    //!           to the held fragment and keep track of the fully connected fragment
    //!
    //! @see @link example_chemistry_fragment_connector.cpp @endlink
    //! @author kothiwsk
    //! @date Oct 10, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API FragmentConnector :
      public util::ObjectInterface
    {
    private:

    //////////
    // data //
    //////////

      //! the current state of fragment
      FragmentComplete                              m_Fragment;

      //! isomorphism between fragment and molecule
      storage::Vector< size_t>                      m_Vertices;

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      FragmentConnector();

      //! @brief Constructor from data members
      //! @param FRAGMENT the base fragment which will be grown by connecting to other fragments
      //! @param VERTICES isomorphism between the base fragment and the molecule which is being assembled
      FragmentConnector
      (
        const FragmentComplete &FRAGMENT,
        const storage::Vector< size_t> &VERTICES
      );

      //! @brief Clone function
      //! @return pointer to new FragmentConnector
      FragmentConnector *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns fragment which this class holds
      //! @return fragment which this class holds
      const FragmentComplete &GetFragment() const
      {
        return m_Fragment;
      }

      //! @brief returns isomorphism of that is being assembled and fragment which this class holds
      //! @return isomorphism of that is being assembled and fragment which this class holds
      const storage::Vector< size_t> &GetVertices() const
      {
        return m_Vertices;
      }

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

    /////////////////
    // operations  //
    /////////////////

      //! @brief merges the incoming fragment to fragment that is already held by the object
      //! @param FRAGMENT the incoming fragment that needs to be merged to the fragment that is already held by the object
      //! @param COMMON_VERTICES atoms that are common between m_Fragment (keys) and incoming fragment(mapped values)
      //! @param ISOMORPHISM isomorphism between incoming fragment and the molecule which is being built
      //! @return true if merge was successful, false otherwise
      bool operator()
      (
        const FragmentComplete &FRAGMENT,
        const storage::Map< size_t, size_t> &COMMON_VERTICES,
        const storage::Vector< size_t> &ISOMORPHISM
      );

      //! @brief merges the incoming fragment to fragment that is already held by the object
      //! @param FRAGMENT the incoming fragment that needs to be added to the fragment that is already held by the object
      //! @param BOND_TYPE the bond that is to be used to connect
      //! @param INDICES_TO_CONNECT atoms of m_Fragment (keys) that need to be connected to atoms of incoming fragment(mapped values)
      //! @param ISOMORPHISM isomorphism between incoming fragment and the molecule which is being built
      //! @return true if merge was successful, false otherwise
      bool operator()
      (
        const FragmentComplete &FRAGMENT,
        const ConfigurationalBondType &BOND_TYPE,
        const storage::Pair< size_t, size_t> &INDICES_TO_CONNECT,
        const storage::Vector< size_t> &ISOMORPHISM
      );

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; // class FragmentConnector

  } // namespace chemistry
} // namespace bcl

#endif // BCL_CHEMISTRY_FRAGMENT_CONNECTOR_H_
