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

// initialize the static initialization fiasco finder, if macro ENABLE_FIASCO_FINDER is defined
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "score/bcl_score_aa_assignments.h"

// includes from bcl - sorted alphabetically
#include "math/bcl_math_z_score.h"
#include "score/bcl_score_aa_assignment_blast_profile.h"
#include "score/bcl_score_aa_assignment_blosum.h"
#include "score/bcl_score_aa_assignment_identity.h"
#include "score/bcl_score_aa_assignment_pam.h"
#include "score/bcl_score_aa_assignment_phat.h"
#include "score/bcl_score_aa_assignment_property.h"
#include "score/bcl_score_aa_assignment_ss_prediction.h"
#include "util/bcl_util_enumerate.hpp"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    AAAssignments::AAAssignments() :
      e_IDENTITY(           AddEnum( "identity"      , util::CloneToShPtr( score::AAAssignmentIdentity()                                                           ))),
      e_PAM100(             AddEnum( "pam100"        , util::CloneToShPtr( score::AAAssignmentPAM( score::AAAssignmentPAM::e_PAM_100)                              ))),
      e_PAM120(             AddEnum( "pam120"        , util::CloneToShPtr( score::AAAssignmentPAM( score::AAAssignmentPAM::e_PAM_120)                              ))),
      e_PAM160(             AddEnum( "pam160"        , util::CloneToShPtr( score::AAAssignmentPAM( score::AAAssignmentPAM::e_PAM_160)                              ))),
      e_PAM250(             AddEnum( "pam250"        , util::CloneToShPtr( score::AAAssignmentPAM( score::AAAssignmentPAM::e_PAM_250)                              ))),
      e_BLOSUM90(           AddEnum( "blosum90"      , util::CloneToShPtr( score::AAAssignmentBLOSUM( score::AAAssignmentBLOSUM::e_BLOSUM_90)                      ))),
      e_BLOSUM80(           AddEnum( "blosum80"      , util::CloneToShPtr( score::AAAssignmentBLOSUM( score::AAAssignmentBLOSUM::e_BLOSUM_80)                      ))),
      e_BLOSUM62(           AddEnum( "blosum62"      , util::CloneToShPtr( score::AAAssignmentBLOSUM( score::AAAssignmentBLOSUM::e_BLOSUM_62)                      ))),
      e_BLOSUM45(           AddEnum( "blosum45"      , util::CloneToShPtr( score::AAAssignmentBLOSUM( score::AAAssignmentBLOSUM::e_BLOSUM_45)                      ))),
      e_PHAT85(             AddEnum( "phat85"        , util::CloneToShPtr( score::AAAssignmentPHAT( score::AAAssignmentPHAT::e_PHAT_85)                            ))),
      e_PHAT80(             AddEnum( "phat80"        , util::CloneToShPtr( score::AAAssignmentPHAT( score::AAAssignmentPHAT::e_PHAT_80)                            ))),
      e_PHAT75(             AddEnum( "phat75"        , util::CloneToShPtr( score::AAAssignmentPHAT( score::AAAssignmentPHAT::e_PHAT_75)                            ))),
      e_PHAT70(             AddEnum( "phat70"        , util::CloneToShPtr( score::AAAssignmentPHAT( score::AAAssignmentPHAT::e_PHAT_70)                            ))),
      e_BLAST(              AddEnum( "blast"         , util::CloneToShPtr( score::AAAssignmentBlastProfile()                                                       ))),
      e_PSIPRED(            AddEnum( "psipred"       , util::CloneToShPtr( score::AAAssignmentSSPrediction( sspred::GetMethods().e_PSIPRED)                        ))),
      e_JUFO(               AddEnum( "jufo"          , util::CloneToShPtr( score::AAAssignmentSSPrediction( sspred::GetMethods().e_JUFO9D)                         ))),
      e_SAM(                AddEnum( "sam"           , util::CloneToShPtr( score::AAAssignmentSSPrediction( sspred::GetMethods().e_SAM)                            ))),
      e_TMHMM(              AddEnum( "tmhmm"         , util::CloneToShPtr( score::AAAssignmentSSPrediction( sspred::GetMethods().e_TMHMM)                          ))),
      e_TMMOD(              AddEnum( "tmmod"         , util::CloneToShPtr( score::AAAssignmentSSPrediction( sspred::GetMethods().e_TMMOD)                          ))),
      e_B2TMPRED(           AddEnum( "b2tmpred"      , util::CloneToShPtr( score::AAAssignmentSSPrediction( sspred::GetMethods().e_B2TMPRED)                       ))),
      e_PROFTMB(            AddEnum( "proftmb"       , util::CloneToShPtr( score::AAAssignmentSSPrediction( sspred::GetMethods().e_PROFTMB)                        ))),
      e_CONPRED(            AddEnum( "conpred"       , util::CloneToShPtr( score::AAAssignmentSSPrediction( sspred::GetMethods().e_CONPRED)                        ))),
      e_STERICAL_PARAMETER( AddEnum( "steric"        , util::CloneToShPtr( score::AAAssignmentProperty( biol::AATypeData::e_StericalParameter)                     ))),
      e_POLARIZABILITY(     AddEnum( "polarizability", util::CloneToShPtr( score::AAAssignmentProperty( biol::AATypeData::e_Polarizability)                        ))),
      e_VOLUME(             AddEnum( "volume"        , util::CloneToShPtr( score::AAAssignmentProperty( biol::AATypeData::e_Volume)                                ))),
      e_HYDROPHOBICITY(     AddEnum( "hydrophobicity", util::CloneToShPtr( score::AAAssignmentProperty( biol::AATypeData::e_Hydrophobicity)                        ))),
      e_ISOELECTRIC_POINT(  AddEnum( "isoelectric"   , util::CloneToShPtr( score::AAAssignmentProperty( biol::AATypeData::e_IsoelectricPoint)                      ))),
      e_TFE_WHITE(          AddEnum( "tfe_white"     , util::CloneToShPtr( score::AAAssignmentProperty( biol::AATypeData::e_TransferFreeEnergyWhimleyWhite)        ))),
      e_TFE_ENGELMAN(       AddEnum( "tfe_engelman"  , util::CloneToShPtr( score::AAAssignmentProperty( biol::AATypeData::e_TransferFreeEnergyEngelmanSeitzGoldman))))
    {
      m_ZScores[ e_IDENTITY          ] = util::ShPtr< math::ZScore>( new math::ZScore(  0.0000, 1.0000));
      m_ZScores[ e_PAM100            ] = util::ShPtr< math::ZScore>( new math::ZScore( -0.1951, 0.3109));
      m_ZScores[ e_PAM120            ] = util::ShPtr< math::ZScore>( new math::ZScore( -0.1581, 0.2803));
      m_ZScores[ e_PAM160            ] = util::ShPtr< math::ZScore>( new math::ZScore( -0.1108, 0.2340));
      m_ZScores[ e_PAM250            ] = util::ShPtr< math::ZScore>( new math::ZScore( -0.0824, 0.2498));
      m_ZScores[ e_BLOSUM90          ] = util::ShPtr< math::ZScore>( new math::ZScore( -0.1698, 0.2586));
      m_ZScores[ e_BLOSUM80          ] = util::ShPtr< math::ZScore>( new math::ZScore( -0.2111, 0.3580));
      m_ZScores[ e_BLOSUM62          ] = util::ShPtr< math::ZScore>( new math::ZScore( -0.0967, 0.2088));
      m_ZScores[ e_BLOSUM45          ] = util::ShPtr< math::ZScore>( new math::ZScore( -0.0821, 0.2273));
      m_ZScores[ e_PHAT85            ] = util::ShPtr< math::ZScore>( new math::ZScore( -0.1431, 0.3114));
      m_ZScores[ e_PHAT80            ] = util::ShPtr< math::ZScore>( new math::ZScore( -0.1676, 0.4248));
      m_ZScores[ e_PHAT75            ] = util::ShPtr< math::ZScore>( new math::ZScore( -0.1469, 0.3990));
      m_ZScores[ e_PHAT70            ] = util::ShPtr< math::ZScore>( new math::ZScore( -0.0895, 0.3552));
      m_ZScores[ e_BLAST             ] = util::ShPtr< math::ZScore>( new math::ZScore( -0.0072, 0.0881));
      m_ZScores[ e_PSIPRED           ] = util::ShPtr< math::ZScore>( new math::ZScore( -0.1431, 0.4728));
      m_ZScores[ e_JUFO              ] = util::ShPtr< math::ZScore>( new math::ZScore( -0.0388, 0.2451));
      m_ZScores[ e_SAM               ] = util::ShPtr< math::ZScore>( new math::ZScore( -0.0056, 0.2076));
      m_ZScores[ e_TMHMM             ] = util::ShPtr< math::ZScore>( new math::ZScore(  0.0000, 0.5000));
      m_ZScores[ e_TMMOD             ] = util::ShPtr< math::ZScore>( new math::ZScore(  0.0000, 0.5000));
      m_ZScores[ e_B2TMPRED          ] = util::ShPtr< math::ZScore>( new math::ZScore(  0.0000, 0.5000));
      m_ZScores[ e_PROFTMB           ] = util::ShPtr< math::ZScore>( new math::ZScore(  0.0000, 0.5000));
      m_ZScores[ e_CONPRED           ] = util::ShPtr< math::ZScore>( new math::ZScore(  0.0000, 0.5000));
      m_ZScores[ e_STERICAL_PARAMETER] = util::ShPtr< math::ZScore>( new math::ZScore( -1.1514, 0.8981));
      m_ZScores[ e_POLARIZABILITY    ] = util::ShPtr< math::ZScore>( new math::ZScore( -0.1061, 0.0814));
      m_ZScores[ e_VOLUME            ] = util::ShPtr< math::ZScore>( new math::ZScore( -1.9938, 1.5660));
      m_ZScores[ e_HYDROPHOBICITY    ] = util::ShPtr< math::ZScore>( new math::ZScore( -1.0737, 0.7871));
      m_ZScores[ e_ISOELECTRIC_POINT ] = util::ShPtr< math::ZScore>( new math::ZScore( -1.6180, 1.8058));
      m_ZScores[ e_TFE_WHITE         ] = util::ShPtr< math::ZScore>( new math::ZScore( -1.8252, 1.4175));
      m_ZScores[ e_TFE_ENGELMAN      ] = util::ShPtr< math::ZScore>( new math::ZScore( -5.3870, 4.7486));

      m_ScoreFile[ e_BLAST   ] = ".ascii";
      m_ScoreFile[ e_PSIPRED ] = ".psipred_ss2";
      m_ScoreFile[ e_JUFO    ] = ".jufo";
      m_ScoreFile[ e_SAM     ] = ".rdb6";
      m_ScoreFile[ e_TMHMM   ] = ".tmhmm";
      m_ScoreFile[ e_TMMOD   ] = ".tmmod";
      m_ScoreFile[ e_B2TMPRED] = ".tmpdb";
      m_ScoreFile[ e_PROFTMB ] = ".proftmb";
      m_ScoreFile[ e_CONPRED ] = ".conpred";
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &AAAssignments::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief z-score associated with that enum
    //! @param ASSIGNMENT_ENUM score
    //! @return z-score to normalize the result of the energy function
    const util::ShPtr< math::ZScore> &AAAssignments::GetZScore( const AAAssignment &ASSIGNMENT_ENUM) const
    {
      // find the zscore for given enum
      const storage::Map< AAAssignment, util::ShPtr< math::ZScore> >::const_iterator itr( m_ZScores.Find( ASSIGNMENT_ENUM));
      if( itr == m_ZScores.End())
      {
        static const util::ShPtr< math::ZScore> s_empty;
        return s_empty;
      }

      // end
      return itr->second;
    }

    //! @brief score file extension, if score requires file to be read
    //! @param ASSIGNMENT_ENUM score
    //! @return the file extension required for the score - empty if non is required
    const std::string &AAAssignments::GetFileExtension( const AAAssignment &ASSIGNMENT_ENUM) const
    {
      // find the file extension
      const storage::Map< AAAssignment, std::string>::const_iterator itr( m_ScoreFile.Find( ASSIGNMENT_ENUM));

      // check that it exists
      if( itr == m_ScoreFile.End())
      {
        static const std::string s_empty;
        return s_empty;
      }

      // return the string
      return itr->second;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief construct on access function for all AAAssignments
    //! @return reference to only instances of AAAssignments
    AAAssignments &GetAAAssignments()
    {
      return AAAssignments::GetEnums();
    }

  } // namespace score

  namespace util
  {

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    template class BCL_API Enumerate< ShPtr< function::BinaryInterface< const biol::AABase, const biol::AABase, double> >, score::AAAssignments>;

  } // namespace util
} // namespace bcl
