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
#include "chemistry/bcl_chemistry_mol_align_by_parts.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_conformation_comparison_property_field_correlation.h"
#include "chemistry/bcl_chemistry_conformation_comparison_psi_flex_field.h"
#include "chemistry/bcl_chemistry_conformation_interface.h"
#include "chemistry/bcl_chemistry_fragment_align_to_scaffold.h"
#include "chemistry/bcl_chemistry_priority_dihedral_angles.h"
#include "chemistry/bcl_chemistry_rotamer_library_file.h"
#include "chemistry/bcl_chemistry_sample_conformations.h"
#include "descriptor/bcl_descriptor_combine.h"
#include "descriptor/bcl_descriptor_constants.h"
#include "io/bcl_io_serialization.h"
#include "util/bcl_util_assert.h"
#include "util/bcl_util_string_functions.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

    // add the interface to the set of known implementations
//    const util::SiPtr< const util::ObjectInterface> MolAlignByParts::s_Instance
//    (
//      util::Enumerated< ConformationComparisonInterface>::AddInstance( new MolAlignByParts())
//    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor
    MolAlignByParts::MolAlignByParts() :
      m_MobileAtoms( storage::Vector< size_t>()),
      m_FixedAtoms( storage::Vector< size_t>()),
      m_TargetMaskedAtoms( storage::Vector< size_t>()),
      m_TargetVisibleAtoms( storage::Vector< size_t>()),
      m_TrialsPerScaffold( 10),
      m_MolAlign( ConformationComparisonPsiField())
    {
    }

    //! virtual copy constructor
    MolAlignByParts *MolAlignByParts::Clone() const
    {
      return new MolAlignByParts( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &MolAlignByParts::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns the name used for this class in an object data label
    //! @return the name used for this class in an object data label
    const std::string &MolAlignByParts::GetAlias() const
    {
      static std::string s_name( "MolAlignByParts");
      return s_name;
    }

    //! @brief returns the atom indices that can be sampled during alignment
    //! @return the mobile atom indices
    const storage::Vector< size_t> &MolAlignByParts::GetMobileAtoms() const
    {
      return m_MobileAtoms;
    }

    //! @brief returns the atom indices that cannot be sampled during alignment
    //! @return the fixed atom indices
    const storage::Vector< size_t> &MolAlignByParts::GetFixedAtoms() const
    {
      return m_FixedAtoms;
    }

    //! @brief returns the atom indices that are hidden from scoring during alignment
    //! @return the masked atom indices
    const storage::Vector< size_t> &MolAlignByParts::GetMaskedAtoms() const
    {
      return m_TargetMaskedAtoms;
    }

    //! @brief returns the atom indices that can scored during alignment
    //! @return the visible atom indices
    const storage::Vector< size_t> &MolAlignByParts::GetVisibleAtoms() const
    {
      return m_TargetVisibleAtoms;
    }

  /////////////////
  //  operations //
  /////////////////

    //! @brief align two small molecule objects and find the property RMSD
    //! @param CONFORMERS_MOL_A - conformers generated from input molecule A
    //! @param CONFORMERS_MOL_B - conformers generated from input molecule B
    storage::Vector< storage::Triplet< FragmentComplete, FragmentEnsemble, double> > MolAlignByParts::Align
    (
      const FragmentComplete &TARGET_MOL,
      const FragmentEnsemble &SCAFFOLD_ENS
    )
    {
      // initialize our return fragment
      FragmentComplete mol( TARGET_MOL);

      // make our mutant scaffold
      FragmentComplete mutant_scaffold( CreateMutantFragment( SCAFFOLD_ENS));

      // setup sample by parts
      InitializeSampleByParts( mol);
      BCL_MessageStd( "Mobile atoms: " + PrintClean( ( m_MobileAtoms)));
      BCL_MessageStd( "Fixed atoms: " + PrintClean( ( m_FixedAtoms)));
      BCL_MessageStd( "Visible atoms: " + PrintClean( ( m_TargetVisibleAtoms)));
      BCL_MessageStd( "Masked atoms: " + PrintClean( ( m_TargetMaskedAtoms)));

      // setup our atom masks
      InitializeAtomMasks( mol, true);
      BCL_MessageStd( "Visible atoms: " + PrintClean( ( m_TargetVisibleAtoms)));
      BCL_MessageStd( "Masked atoms: " + PrintClean( ( m_TargetMaskedAtoms)));

      // remember the original params
      storage::Vector< size_t> mobile_atoms( m_MobileAtoms);
      storage::Vector< size_t> fixed_atoms( m_FixedAtoms);
      storage::Vector< size_t> target_masked_atoms( m_TargetMaskedAtoms);
      storage::Vector< size_t> target_visible_atoms( m_TargetVisibleAtoms);

      // now what?
      // start with original use-case, then generalize
      // a. make confs
      // b. find best scoring conf across a scaffold (start with 1)
      // c. get mutually matched atoms, invert, find atoms without matches
      // d. set mobile atoms and visible atoms to those that do NOT have matches in the scaffold
      // e. loop back at (a)

      // initialize conformer generator
      static RotamerLibraryFile rotlib;
      static SampleConformations sample_confs
      (
        rotlib,                     // rotamer library file
        "SymmetryRMSD",             // conformation comparer type
        0.25,                       // conformational comparer tolerance
        500,                        // number of conformations
        2500,                       // number of iterations
        false,                      // no change chirality
        0.0,                        // random dihedral change weight
        false,                      // generate 3d?
        0.1,                        // clash tolerance
        true                        // cluster?
      );
      bool local( false);
      if( local)
      {
        // only sample bond angles/lengths
        sample_confs.SetSamplingPreferences( false, false, true, false);
      }

        // alignment properties and weights
        util::ObjectDataLabel property_labels
        (
          "("
          "Combine(Define(ChargeNegative=Multiply(Less(lhs=Atom_SigmaCharge,rhs=Constant(0)),Atom_SigmaCharge)),"
          "Define(ChargePositive=Multiply(Greater(lhs=Atom_SigmaCharge,rhs=Constant(0)),Atom_SigmaCharge)),"
          "Define(HbondAcceptorsStrict=Multiply(Atom_HbondAcceptors,Not(Atom_HbondDonors),LessEqual(lhs=BondTypeCount,rhs=Constant(2)))),"
          "Define(ENegOffset=Subtract(lhs=Atom_ElectroNegativity,rhs=Constant(2.5))),"
          "Define(IsENeg=Greater(lhs=ENegOffset,rhs=Constant(0))),"
          "Define(PolarTernary=Subtract(lhs=Multiply(Add(HbondAcceptorsStrict, Atom_HbondDonors), Constant(2)), rhs=Constant(1))),"
          "Define(Atom_Hydrophobic=Not(DescriptorSum(Abs(Partial(indices(1,2,5,8),2DASign(property=PolarTernary, steps=3)))))),"
          "Atom_SigmaCharge,"
          "HbondAcceptorsStrict,"
          "Atom_HbondDonors,"
          "Atom_Polarizability,"
          "Atom_AromaticityAxes,"
          "Multiply(IsENeg,ENegOffset),"
          "Atom_Hydrophobic,"
          "Atom_VDWVolume)"
          ")"
        );
      descriptor::Combine< AtomConformationalInterface, float> properties( PrepareProperties( property_labels));
      storage::Vector< float> property_weights
      (
        util::SplitStringToNumerical< float>( util::TrimString( "5,3.5,7,1.71,2.41,2.41,2.41,3.0,2.0,0.714"), " \t\n\r,")
      );

      // increment max atom distance over course of alignment
      size_t n_trials( 10);
      float dmax_increment( ( 1.15 - 0.70) / ( n_trials - 1));
      storage::Vector< float> max_atom_distances( n_trials);
      for( size_t i( 0), sz( max_atom_distances.GetSize()); i < sz; ++i)
      {
        max_atom_distances( i) = 0.70 + ( dmax_increment * i);
      }

      // initialize our alignment scorer with the properties and weights and whatnot
      ConformationComparisonPropertyFieldCorrelation alignment_scorer
      (
        properties,
        property_weights,
        0.70,
        false,
        1.0e-2,
        0.6,
        2.0,
        1.0e-2
      );
      alignment_scorer.SetPropertyWeights( property_weights);
      //BCL_Debug( alignment_scorer.GetMaxAtomDistance());

      // make movie
//      io::OFStream debug_movie;
//      io::File::MustOpenOFStream(debug_movie, "molalign_by_parts_test.movie.sdf");

      // align to each scaffold
      BCL_MessageStd( "SampleByParts with the following atom indices: " + util::Format()( mol.GetMDLProperty( "SampleByParts")));
      storage::Triplet< FragmentComplete, FragmentComplete, double> best_pose;
      best_pose.Third() = double( 999.9);
//      for
//      (
//          auto scaffold_itr( SCAFFOLD_ENS.Begin()), scaffold_itr_end( SCAFFOLD_ENS.End());
//          scaffold_itr != scaffold_itr_end;
//          ++scaffold_itr
//      )
//      {
        // multiple rounds for each scaffold
        for( size_t round( 0); round < n_trials; ++round)
        {
          // set max atom distance
          alignment_scorer.SetMaxAtomDistance( max_atom_distances( round));
          //BCL_Debug( alignment_scorer.GetMaxAtomDistance());
          //BCL_Debug( max_atom_distances( round));

          //! molalign object if alignment is to allow movement of more than internal DOFs
          ConformationComparisonPsiField m_MolAlign;

          // make conformers
          FragmentEnsemble confs( sample_confs( mol).First());
//          io::OFStream debug_confs;
//          io::File::MustOpenOFStream(debug_confs, "molalign_by_parts_test.confs.sdf");
//          confs.WriteMDL( debug_confs);
//          io::File::CloseClearFStream( debug_confs);

          // need to re-align the fixed atoms to their original pose
          FragmentAlignToScaffold ats
          (
            ConformationGraphConverter::e_AtomType,
            ConfigurationalBondTypeData::e_BondOrderAmideWithIsometryOrAromaticWithRingness,
            size_t( 3)
          );

          // score conformers against scaffold
          for
          (
              auto conf_itr( confs.Begin()), conf_itr_end( confs.End());
              conf_itr != conf_itr_end;
              ++conf_itr
          )
          {
            // align back to start coordinates
            bool self( true);
            if( m_FixedAtoms.GetSize())
            {
              self ?
                  ats.AlignToScaffold( *conf_itr, mol, fixed_atoms, fixed_atoms) :
                  ats.AlignToScaffold( *conf_itr, mutant_scaffold);
            }

            // score conformer
            auto conf_aln_score
            (
              alignment_scorer.PropertyDistanceScorer
              (
                *conf_itr,
                mutant_scaffold,
                1.0,
//                max_atom_distances( round),
                m_TargetMaskedAtoms
              )
            );

            // save if improved
            //BCL_Debug( conf_aln_score.First());
            //BCL_Debug( conf_aln_score.Second());
            if( conf_aln_score.First() < best_pose.Third())
            {
              best_pose.First() = *conf_itr;
              best_pose.Second() = mutant_scaffold;
              best_pose.Third() = conf_aln_score.First();
              //BCL_Debug( best_pose.Third());
              //BCL_Debug( conf_aln_score.Second());
//              best_pose.First().WriteMDL( debug_movie);
            }
          }
          // end of round; get matched atoms
          storage::Vector< size_t> scaffold_atom_indices;
          for( size_t s_i( 0), s_sz( best_pose.Second().GetSize()); s_i < s_sz; ++s_i)
          {
            scaffold_atom_indices.PushBack( s_i);
          }
          auto mutually_matched_atoms( alignment_scorer.GetAlignedAtoms( best_pose.First(), best_pose.Second(), m_TargetVisibleAtoms, scaffold_atom_indices));
          BCL_MessageStd( "Target mutually matched atoms: " + PrintClean( ( mutually_matched_atoms.First())));
          BCL_MessageStd( "Scaffold mutually matched atoms: " + PrintClean( ( mutually_matched_atoms.Second())));

          // add the matched atoms to the fixed atom list
          for( size_t i( 0), sz( mutually_matched_atoms.First().GetSize()); i < sz; ++i)
          {
            m_FixedAtoms.PushBack( mutually_matched_atoms.First()( i));
          }

          // reset the mobile atoms; keep initial atoms rigid
          m_MobileAtoms.Reset();
          InitializeSampleByParts( mol);
          InitializeAtomMasks(mol, true);
          BCL_MessageStd( "Mobile atoms: " + PrintClean( ( m_MobileAtoms)));
          BCL_MessageStd( "Fixed atoms: " + PrintClean( ( m_FixedAtoms)));
          BCL_MessageStd( "Visible atoms: " + PrintClean( ( m_TargetVisibleAtoms)));
          BCL_MessageStd( "Masked atoms: " + PrintClean( ( m_TargetMaskedAtoms)));

          // logical termination
          if( !m_MobileAtoms.GetSize())
          {
            break;
          }
        }
//        // reset stuff
//        m_MobileAtoms = mobile_atoms;
//        m_FixedAtoms = fixed_atoms;
//        m_TargetMaskedAtoms = target_masked_atoms;
//        m_TargetVisibleAtoms = target_visible_atoms;
//      }

//      io::File::CloseClearFStream( debug_movie);
      //BCL_Debug( best_pose.Third());
      io::OFStream debug_out;
      io::File::MustOpenOFStream( debug_out, "molalign_by_parts_test.sdf", std::ios::app);
      best_pose.First().WriteMDL( debug_out);
      io::File::CloseClearFStream( debug_out);

      return storage::Vector< storage::Triplet< FragmentComplete, FragmentEnsemble, double> >();
    }

    //! @brief set atom properties and remove size zero properties
    //! @param PROPERTIES descriptor label to use
    descriptor::Combine< AtomConformationalInterface, float> MolAlignByParts::PrepareProperties( const util::ObjectDataLabel &PROPERTIES) const
    {
      descriptor::Combine< AtomConformationalInterface, float> properties, read_check;

      // Properties must be molecule-wide; assert they are valid
      read_check.SetDimension( 0);
      read_check.AssertRead( PROPERTIES);

      // Fetch the aliases for each property
      for
      (
        descriptor::Combine< AtomConformationalInterface, float>::const_iterator
        itr_prop( read_check.Begin()), itr_prop_end( read_check.End());
        itr_prop != itr_prop_end;
        ++itr_prop
      )
      {
        // Do not include 0-size descriptors, e.g. Define, ForEach, etc...
        if( ( *itr_prop)->GetSizeOfFeatures())
        {
          BCL_MessageVrb( "Preparing properties: " + itr_prop->GetString());
          properties.PushBack( *itr_prop);
        }
      }

      return properties;
    } // PrepareProperties

    //! @brief initialize the sample by parts atoms from member data
    void MolAlignByParts::InitializeSampleByParts( FragmentComplete &MOL)
    {
      storage::Set< size_t> unique_mobile_atoms;
      storage::Vector< size_t> pruned_mobile_atoms;

      // give preference to specified mobile atoms
      if( m_MobileAtoms.GetSize())
      {
        // get unique atoms
        for
        (
            auto itr( m_MobileAtoms.Begin()), itr_end( m_MobileAtoms.End());
            itr != itr_end;
            ++itr
        )
        {
          unique_mobile_atoms.Insert( *itr);
        }
      }
      // if no mobile atoms specified, get the inverse of any specified fixed atoms
      else if( m_FixedAtoms.GetSize())
      {
        // collect mobile atom indices
        storage::Vector< size_t> mobile_atoms;
        for( size_t i( 0), sz( MOL.GetSize()); i < sz; ++i)
        {
          if( !( m_FixedAtoms.Find( i) < m_FixedAtoms.GetSize()))
          {
            mobile_atoms.PushBack( i);
          }
        }

        // get unique atoms
        for
        (
            auto itr( mobile_atoms.Begin()), itr_end( mobile_atoms.End());
            itr != itr_end;
            ++itr
        )
        {
          unique_mobile_atoms.Insert( *itr);
        }

      }
      // set all atoms to mobile by default
      else
      {
        for( size_t i( 0), sz( MOL.GetSize()); i < sz; ++i)
        {
          unique_mobile_atoms.Insert( i);
        }
      }

      // set the SampleByParts property for our molecule
      pruned_mobile_atoms = storage::Vector< size_t>( unique_mobile_atoms.Begin(), unique_mobile_atoms.End());
      MOL.GetStoredPropertiesNonConst().SetMDLProperty( "SampleByParts", pruned_mobile_atoms);
      BCL_MessageStd( "SampleByParts with the following atom indices: " + util::Format()( MOL.GetMDLProperty( "SampleByParts")));

      // update member data
      m_MobileAtoms = pruned_mobile_atoms;
      m_FixedAtoms.Reset();
      for( size_t i( 0), sz( MOL.GetSize()); i < sz; ++i)
      {
        if( !( m_MobileAtoms.Find( i) < m_MobileAtoms.GetSize()))
        {
          m_FixedAtoms.PushBack( i);
        }
      }
    }

    //! @brief initialize masking from member data
    //! @param MOL - target molecule
    //! @param MASK_FIXED_ATOMS - only allow scoring on mobile atoms; mask the fixed atoms
    void MolAlignByParts::InitializeAtomMasks
    (
      const FragmentComplete &MOL,
      const bool &MASK_FIXED_ATOMS
    )
    {
      // run InitializeSampleByParts first, then just set it so that scoring only occurs on mobile atoms
      if( MASK_FIXED_ATOMS && m_MobileAtoms.GetSize() && m_FixedAtoms.GetSize())
      {
        m_TargetVisibleAtoms = m_MobileAtoms;
        m_TargetMaskedAtoms = m_FixedAtoms;
      }
      // reference masked atoms
      else if( m_TargetMaskedAtoms.GetSize())
      {
        m_TargetVisibleAtoms.Reset();
        for( size_t i( 0), sz( MOL.GetSize()); i < sz; ++i)
        {
          // do not include the atoms we said we wanted to exclude
          bool exclude( false);
          for( size_t e_i( 0), e_sz( m_TargetMaskedAtoms.GetSize()); e_i < e_sz; ++e_i)
          {
            if( i == m_TargetMaskedAtoms( e_i))
            {
              exclude = true;
              break;
            }
          }
          // if we did not find it in the exclusion indices, add it
          if( !exclude)
          {
            m_TargetVisibleAtoms.PushBack( i);
          }
        }
      }
      // reference visible atoms
      else if( m_TargetVisibleAtoms.GetSize())
      {
        m_TargetMaskedAtoms.Reset();
        for( size_t i( 0), sz( MOL.GetSize()); i < sz; ++i)
        {
          // do not include the atoms we said we wanted to exclude
          bool exclude( false);
          for( size_t e_i( 0), e_sz( m_TargetVisibleAtoms.GetSize()); e_i < e_sz; ++e_i)
          {
            if( i == m_TargetVisibleAtoms( e_i))
            {
              exclude = true;
              break;
            }
          }
          // if we did not find it in the exclusion indices, add it
          if( !exclude)
          {
            m_TargetMaskedAtoms.PushBack( i);
          }
        }
      }
      // allow scoring on all atoms
      else
      {
        for( size_t i( 0), sz( MOL.GetSize()); i < sz; ++i)
        {
          m_TargetVisibleAtoms.PushBack( i);
        }
      }

      // finalize by filtering indices for uniqueness
      if( m_TargetVisibleAtoms.GetSize())
      {
        storage::Set< size_t> unique_atoms;
        for( size_t i( 0), sz( m_TargetVisibleAtoms.GetSize()); i < sz; ++i)
        {
          unique_atoms.Insert( m_TargetVisibleAtoms( i));
        }
        m_TargetVisibleAtoms = storage::Vector< size_t>( unique_atoms.Begin(), unique_atoms.End());
      }
      if( m_TargetMaskedAtoms.GetSize())
      {
        storage::Set< size_t> unique_atoms;
        for( size_t i( 0), sz( m_TargetMaskedAtoms.GetSize()); i < sz; ++i)
        {
          unique_atoms.Insert( m_TargetMaskedAtoms( i));
        }
        m_TargetMaskedAtoms = storage::Vector< size_t>( unique_atoms.Begin(), unique_atoms.End());
      }
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief set the atom indices that can be sampled during alignment
    void MolAlignByParts::SetMobileAtoms( const storage::Vector< size_t> &MOBILE_ATOM_INDICES)
    {
      m_MobileAtoms = MOBILE_ATOM_INDICES;
    }

    //! @brief set the atom indices that are fixed during alignment
    void MolAlignByParts::SetFixedAtoms( const storage::Vector< size_t> &FIXED_ATOM_INDICES)
    {
      m_FixedAtoms = FIXED_ATOM_INDICES;
    }

    //! @brief set the atom indices that are hidden from scoring during alignment
    void MolAlignByParts::SetMaskedAtoms( const storage::Vector< size_t> &MASKED_ATOM_INDICES)
    {
      m_TargetMaskedAtoms = MASKED_ATOM_INDICES;
    }

    //! @brief set the atom indices that can be scored during alignment
    void MolAlignByParts::SetVisibleAtoms( const storage::Vector< size_t> &VISIBLE_ATOM_INDICES)
    {
      m_TargetVisibleAtoms = VISIBLE_ATOM_INDICES;
    }

    //! @brief combine a bunch of fragments without destroying internal connectivity
    FragmentComplete MolAlignByParts::CreateMutantFragment( const FragmentEnsemble &ENSEMBLE) const
    {
      // Initialize sdf info with data from first molecule in ensemble
      storage::Vector< sdf::AtomInfo> new_atominfo( ENSEMBLE.GetMolecules().FirstElement().GetAtomInfo());
      storage::Vector< sdf::BondInfo> new_bondinfo( ENSEMBLE.GetMolecules().FirstElement().GetBondInfo());

      // iterate over remainder of molecules
      size_t mol_index( 0);
      for
      (
          auto mol_itr( ENSEMBLE.Begin()), mol_itr_end( ENSEMBLE.End());
          mol_itr != mol_itr_end;
          ++mol_itr, ++mol_index
      )
      {
        // skip first molecule in ensemble; already added
        if( mol_index == size_t( 0))
        {
          continue;
        }

        // keep track of molecule size
        size_t current_mol_size( new_atominfo.GetSize());

        // Add atominfo from this molecule to new atominfo
        storage::Vector< sdf::AtomInfo> atominfo( mol_itr->GetAtomInfo());
        for( size_t atominfo_index( 0), atominfo_sz( atominfo.GetSize()); atominfo_index < atominfo_sz; ++atominfo_index)
        {
          new_atominfo.PushBack( atominfo( atominfo_index));
        }

        // Add bondinfo from this molecule to new bond info
        storage::Vector< sdf::BondInfo> bondinfo( mol_itr->GetBondInfo());
        for( size_t bondinfo_index( 0), bondinfo_sz( bondinfo.GetSize()); bondinfo_index < bondinfo_sz; ++bondinfo_index)
        {
          new_bondinfo.PushBack
          (
            sdf::BondInfo
            (
              bondinfo( bondinfo_index).GetAtomIndexLow() + current_mol_size,
              bondinfo( bondinfo_index).GetAtomIndexHigh() + current_mol_size,
              bondinfo( bondinfo_index).GetConfigurationalBondType()
            )
          );
        }
      }

      // return the new molecule
      AtomVector< AtomComplete> new_vector( new_atominfo, new_bondinfo);
      return FragmentComplete( new_vector, "");
    }

    //! @brief write a BCL vector in a horizontal string format for easier reading
    std::string MolAlignByParts::PrintClean( const storage::Vector< size_t> &INPUT) const
    {
      std::stringstream ostream;
      for
      (
          auto itr( INPUT.Begin()), itr_end( INPUT.End());
          itr != itr_end;
          ++itr
      )
      {
        ostream << util::Format()( *itr);
        ostream << " ";
      }
      return ostream.str();
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &MolAlignByParts::Read( std::istream &ISTREAM)
    {
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &MolAlignByParts::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      return OSTREAM;
    }

  } // namespace chemistry
} // namespace bcl
