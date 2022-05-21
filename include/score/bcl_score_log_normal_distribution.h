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

#ifndef BCL_SCORE_LOG_NORMAL_DISTRIBUTION_H_
#define BCL_SCORE_LOG_NORMAL_DISTRIBUTION_H_

// include the namespace header
#include "bcl_score.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_score_restraint_nmr_distance_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class LogNormalDistribution
    //! @brief Provides another way to score distance restraints.
    //! @details ((log( obs_dist / calc_dist))^2)-1
    //!        Nilges Et. Al. 2005 JACS "Modeling Errors in NOE Data with a Log-normal Distribution Improves the Quality
    //!        of NMR Structures"
    //!
    //! @see @link example_score_log_normal_distribution.cpp @endlink
    //! @author akinlr, alexanns
    //! @date Jul 9, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API LogNormalDistribution :
      public RestraintNMRDistanceInterface
    {

    private:

    //////////
    // data //
    //////////

      static const double s_WellDepth; //!< where the minima of the graph should fall

      double m_KVariable;             //!< defines harshness of scoring function

      //! scheme to be used
      std::string m_Scheme;

      //! score for a restraint with residues/atoms not found in the protein model
      static const double s_DefaultScore;

      //! effective distance per bond
      static const double s_EffectiveDistancePerBond;

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
      LogNormalDistribution();

      //! @brief parameter constructor
      //! @param KVARIABLE
      explicit LogNormalDistribution( const double &KVARIABLE);

      //! @brief Clone function
      //! @return pointer to new LogNormalDistribution
      virtual LogNormalDistribution *Clone() const;

      //! @brief virtual destructor
      virtual ~LogNormalDistribution();

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief returns default scheme
      //! @return default scheme
      static const std::string &GetDefaultScheme();

      //! @brief returns scheme being used
      //! @return scheme being used
      const std::string &GetScheme() const;

    ///////////////
    // operators //
    ///////////////

      //! @brief () operator scores protein model
      //! @param RESTRAINT restraint to be scored
      //! @return score
      double operator()( const restraint::AtomDistanceAssignment &RESTRAINT) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    }; // class LogNormalDistribution

  } // namespace score
} // namespace bcl

#endif // BCL_SCORE_LOG_NORMAL_DISTRIBUTION_H_
