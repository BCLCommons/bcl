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

#ifndef BCL_CHEMISTRY_ATOM_ENVIRONMENT_H_
#define BCL_CHEMISTRY_ATOM_ENVIRONMENT_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_atom_conformational_interface.h"
#include "storage/bcl_storage_map.h"
#include "util/bcl_util_sh_ptr.h"
#include "util/bcl_util_si_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AtomEnvironment
    //! @brief stores the atoms connected to one substituent of the 'atom of interest'. This class is used for
    //! describing the environment of each 'atom of interest' by generating a cone. Each substituent is encoded in
    //! its own cone.
    //!
    //! @see @link example_chemistry_atom_environment.cpp @endlink
    //! @author mueller, mendenjl, brownbp1
    //! @date 08/24/2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API AtomEnvironment :
      public util::ObjectInterface
    {

    private:

    //////////
    // data //
    //////////

      //! the atom of interest
      util::SiPtr< const AtomConformationalInterface>   m_AtomOfInterest;

      //! the actual atom environment
      storage::Map< util::SiPtr< const AtomConformationalInterface>, storage::Pair< size_t, double> > m_EnvironmentAtoms;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! default constructor
      AtomEnvironment();

      //! construct from number of spheres
      AtomEnvironment
      (
        const AtomConformationalInterface &ATOM_OF_INTEREST,
        const size_t NUMBER_SPHERES
      );

      //! virtual copy constructor
      AtomEnvironment *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief get the map that describes the environment of the atom of interest
      //! @return the map that describes the environment of the atom of interest
      const storage::Map< util::SiPtr< const AtomConformationalInterface>, storage::Pair< size_t, double> > &GetEnvironmentAtoms() const
      {
        return m_EnvironmentAtoms;
      }

      //! @brief compute atom environment two bonds from reference atom
      //! @param A atom of interest
      //! @return return string indicating atom environment
      static std::string MakeAtomEnvironmentStringTwo( const AtomConformationalInterface &A);

      //! @brief compute atom environment three bonds from reference atom
      //! @param A atom of interest
      //! @return return string indicating atom environment
      static std::string MakeAtomEnvironmentStringThree( const AtomConformationalInterface &A);

      //! @brief compute atom environment four bonds from reference atom
      //! @param A atom of interest
      //! @return return string indicating atom environment
      static std::string MakeAtomEnvironmentStringFour( const AtomConformationalInterface &A);

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! write AtomEnvironment to std::ostream
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

      //! read AtomEnvironment from io::IFStream
      std::istream &Read( std::istream &ISTREAM);

    }; // class AtomEnvironment

  } // namespace chemistry
} // namespace bcl

#endif // BCL_CHEMISTRY_ATOM_ENVIRONMENT_H_
