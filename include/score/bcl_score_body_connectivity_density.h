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

#ifndef BCL_SCORE_BODY_CONNECTIVITY_DENSITY_H_
#define BCL_SCORE_BODY_CONNECTIVITY_DENSITY_H_

// include the namespace header
#include "bcl_score.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"
#include "density/bcl_density.fwd.hh"
#include "restraint/bcl_restraint.fwd.hh"
#include "storage/bcl_storage.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_score_protein_model.h"
#include "density/bcl_density_connectivity.h"
#include "restraint/bcl_restraint_body.h"
#include "storage/bcl_storage_map.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class BodyConnectivityDensity
    //! @brief is for scoring the strength of density that is between two ends of restraint
    //!        Bodies. This is based on an experimental density map.
    //!        If there is strong density between the two bodies, it is likely that these two bodies are
    //!        connected in sequence. Thus a protein model connecting these two bodies (density rods) should get a
    //!        better score.
    //!
    //! @see @link example_score_body_connectivity_density.cpp @endlink
    //! @author alexanns, linders
    //! @date February 16, 2009
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API BodyConnectivityDensity :
      public ProteinModel
    {

    private:

    //////////
    // data //
    //////////

      //! storage::Map "m_Scores" contains the score for connecting any of the restraint::Bodies
      storage::Map< density::Connectivity, double, density::Connectivity::LessThan> m_Scores;

      //! ShPtr< restraint::Body> "m_BodyRestraint" contains the information about bodies between which density scores
      //! are calculated
      util::ShPtr< restraint::Body> m_BodyRestraint;

      //! scheme to be used in outputting schemes
      std::string m_Scheme;

    //////////
    // data //
    //////////

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

      //! @brief returns default scheme
      //! @return default scheme
      static const std::string &GetDefaultScheme();

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      BodyConnectivityDensity();

      //! @brief constructor taking a storage::Map and  ShPtr< restraint::Body>
      //! @param CONNECTIVITY_SCORES the storage::Map which will be "m_Scores"
      //! @param BODY_RESTRAINT util::ShPtr< restraint::Body> which will be "m_BodyRestraint"
      //! @param SCHEME
      BodyConnectivityDensity
      (
        const storage::Map< density::Connectivity, double, density::Connectivity::LessThan> &CONNECTIVITY_SCORES,
        const util::ShPtr< restraint::Body> &BODY_RESTRAINT,
        const std::string &SCHEME = GetDefaultScheme()
      );

      //! @brief constructor taking a restraint::Body and a Density Map
      //!        This constructor calculates the scores that go into "m_Scores" and fills "m_Scores"
      //! @param BODY_RESTRAINT contains the information about which bodies densities will be calculated between and
      //!        which will be "m_BodyRestraint"
      //! @param DENSITY_MAP experimental density map which is used to calculate the minimal density between each of
      //!        the bodies in "BODY_RESTRAINT"
      //! @param SCHEME
      BodyConnectivityDensity
      (
        const util::ShPtr< restraint::Body> &BODY_RESTRAINT,
        const util::ShPtr< density::Map> &DENSITY_MAP,
        const std::string &SCHEME = GetDefaultScheme()
      );

      //! @brief copy constructor
      BodyConnectivityDensity *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief gets the scheme for this function class
      //! @return the scheme for this function class
      const std::string &GetScheme() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief GetScores gives a const reference to "m_Scores"
      //! @return a const reference to "m_Scores"
      const storage::Map< density::Connectivity, double, density::Connectivity::LessThan> &GetScores() const;

      //! @brief GetBodyRestraint returns a const reference to "m_BodyRestraint"
      //! @return a const reference to "m_BodyRestraint"
      const util::ShPtr< restraint::Body> &GetBodyRestraint() const;

    ///////////////
    // operators //
    ///////////////

      //! @brief operator() which takes an Assignment and returns the density connectivity score.
      //!        The assignment contains all the body restraints and the corresponding protein model bodies.
      //!        For each pair of restraint bodies which are related by sequence (i.e. would be connected by a loop)
      //!        according to the protein model bodies which fill them,
      //!        the two points between which the density needs to be calculated are determined. The score for the
      //!        density between the two points is then determined.
      //! @param PROTEIN_MODEL the ProteinModel which will be scored for agreement with the connectivity
      //! @return return a double which is the score of the agreement of the ProteinModel with the connectivity
      double operator()
      (
        const assemble::ProteinModel &PROTEIN_MODEL
      ) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read distance from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write distance to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT indentation
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    //////////////////////
    // helper functions //
    //////////////////////

    protected:

      //! @brief InitializeScores initializes "m_Scores" upon construction
      //! @param BODY_RESTRAINT contains the information about which bodies densities will be calculated between and
      //!        which will be "m_BodyRestraint"
      //! @param DENSITY_MAP experimental density map which is used to calculate the minimal density between each of
      //!        the bodies in "BODY_RESTRAINT"
      //! @return storage::Map which will be "m_Scores"
      storage::Map< density::Connectivity, double, density::Connectivity::LessThan>
      InitializeScores
      (
        const util::ShPtr< restraint::Body> &BODY_RESTRAINT,
        const util::ShPtr< density::Map> &DENSITY_MAP
      ) const;

      //! @brief ConvertIntensitiesToScores takes the raw intensities determined between restraint::Bodies and
      //!        converts them into a Z-score
      //! @param CONNECTIVITIES a storage::List which has all the DensityConnectivity objects to be converted to scores.
      //! @return returns a storage::Map which has the Z-scores and will be "m_Scores"
      storage::Map< density::Connectivity, double, density::Connectivity::LessThan>
      ConvertIntensitiesToScores
      (
        const storage::List< density::Connectivity> &CONNECTIVITIES
      ) const;

      //! @brief GetRestraintBodiesandOrientations determines the two restraint bodies (i.e. density rods) occupied by
      //!        two neighboring SSEs in the protein model and at which ends of the restraint bodies the protein model
      //!        connects them (in the form of two bools)
      //! @param SSES are the two SSEs for which the two points of connection are desired
      //! @return returns a storage::Pair of vector of two bodies and vector of two bools
      storage::Pair< storage::VectorND< 2, const util::ShPtr< coord::GeometryInterface> >, storage::VectorND< 2, bool> >
      GetRestraintBodiesandOrientations( const storage::VectorND< 2, util::SiPtr< const assemble::SSE> > &SSES) const;

      //! @brief DetermineSSEOrientationInBodyRestraint is used to match a coordinate with one of the end points of
      //!        a body by checking which end the coordinate is closer to. This is used to determine the which
      //!        direction an SSE is pointing in a restraining body, for example.
      //! @param SSE_POINT is the coordinate for which its position in the restraining body is desired
      //! @param BODY the body for which the placement position of SSE_POINT is desired
      //! @return returns a bool to indicate to which of the ends of "BODY" "SSE_POINT" is closest (true for beginning
      //!         of BODY, false for end of BODY)
      bool DetermineSSEOrientationInBodyRestraint( const linal::Vector3D &SSE_POINT, const coord::GeometryInterface &BODY) const;

      //! @brief GetScore returns the calculated Z-score connectivity score for connecting two bodies (at specific ends)
      //! @param CONNECTIVITY_INFORMATION a storage::Pair of vector of two bodies and vector of two bools (specifying
      //!        which two bodies are to be connected and at which ends of the bodies)
      //! @return returns a double which is the Z-score associated with connecting CONNECTIVITY_INFORMATION
      double CalculateScore
      (
        const storage::Pair< storage::VectorND< 2, const util::ShPtr< coord::GeometryInterface> >, storage::VectorND< 2, bool> > &
        CONNECTIVITY_INFORMATION
      ) const;

    }; // class BodyConnectivityDensity

  } // namespace score
} // namespace bcl

#endif //BCL_SCORE_BODY_CONNECTIVITY_DENSITY_H_
