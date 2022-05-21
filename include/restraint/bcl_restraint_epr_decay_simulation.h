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

#ifndef BCL_RESTRAINT_EPR_DECAY_SIMULATION_H_
#define BCL_RESTRAINT_EPR_DECAY_SIMULATION_H_

// include the namespace header
#include "bcl_restraint.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "storage/bcl_storage_vector.h"
#include "util/bcl_util_serializable_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class EPRDecaySimulation
    //! @brief Interface for optimization methods
    //!
    //! @see @link example_restraint_epr_decay_simulation.cpp @endlink
    //! @author fischea
    //! @date Nov 3, 2016
    //!
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////
    class BCL_API EPRDecaySimulation :
      public util::SerializableInterface
    {

    //////////////
    // typedefs //
    //////////////

    public:

      //! representation of one spin-labeling site consisting of chain ID and sequence ID
      typedef storage::Pair< char, int> SLSite;

      //! representation of one spin-labeling pair
      typedef storage::Pair< SLSite, SLSite> SLPair;

    //////////
    // data //
    //////////

      //! single instance of this class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    private:

      //! list of pairs of spin-labeling sites
      storage::Vector< SLPair> m_SpinLabelingPairs;

    ////////////////////////////////////
    //  construction and destruction  //
    ////////////////////////////////////

    public:

      //! @brief default constructor
      EPRDecaySimulation();

      //! @brief construct from list of spin-labeling pairs
      //! @param SL_PAIRS list of spin-labeling pairs
      EPRDecaySimulation( const storage::Vector< SLPair> &SL_PAIRS);

      //! @brief virtual copy constructor
      //! @return pointer to a new EPRDecaySimulation
      EPRDecaySimulation *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns the name of this class
      //! @return the name of this class
      const std::string &GetClassIdentifier() const;

      //! @brief get the name of the object when used in a dynamic context
      //! @return the name of the object when used in a dynamic context
      const std::string &GetAlias() const;

      //! @brief returns parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief calculate decay for the provided distance and time
      //! @param DISTANCE spin-spin distance
      //! @param TIME
      //! @return decay
      static double CalculateDecay( double DISTANCE, double TIME);

    ///////////////
    // operators //
    ///////////////

    }; // class EPRDecaySimulation

  } // namespace restraint
} // namespace bcl

#endif // BCL_RESTRAINT_EPR_DECAY_SIMULATION_H_
