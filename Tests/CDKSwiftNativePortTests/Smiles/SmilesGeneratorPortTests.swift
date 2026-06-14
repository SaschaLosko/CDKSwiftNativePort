import XCTest
@testable import CDKSwiftNativePort

final class SmilesGeneratorPortTests: XCTestCase {
    private let parser = CDKSmilesParserFactory.shared.newSmilesParser(flavor: .cdkDefault)

    func testGeneratesSmilesForEthanol() throws {
        let molecule = try parser.parseSmiles("CCO")
        let generator = CDKSmilesGeneratorFactory.shared.newSmilesGenerator(flavor: [.useAromaticSymbols, .strict])

        let smiles = generator.create(molecule)
        XCTAssertEqual(smiles, "CCO")
    }

    func testGeneratesIsomericSmilesWithChirality() throws {
        let molecule = try parser.parseSmiles("C[C@H](N)O")
        let generator = CDKSmilesGeneratorFactory.shared.newSmilesGenerator(
            flavor: [.useAromaticSymbols, .isomeric, .strict]
        )

        let isomeric = generator.create(molecule)
        XCTAssertTrue(isomeric.contains("@"))

        let reparsed = try parser.parseSmiles(isomeric)
        XCTAssertEqual(reparsed.atomCount, molecule.atomCount)
        XCTAssertEqual(reparsed.bondCount, molecule.bondCount)
    }

    func testGeneratesRingClosureNotation() throws {
        let molecule = try parser.parseSmiles("c1ccccc1")
        let generator = CDKSmilesGeneratorFactory.shared.newSmilesGenerator(
            flavor: [.useAromaticSymbols, .isomeric, .strict]
        )

        let smiles = generator.create(molecule)
        XCTAssertTrue(smiles.contains("1"))

        let reparsed = try parser.parseSmiles(smiles)
        XCTAssertEqual(reparsed.atomCount, 6)
        XCTAssertEqual(reparsed.bondCount, 6)
    }

    func testGeneratesDisconnectedComponents() throws {
        let molecule = try parser.parseSmiles("[Na+].[OH-]")
        let generator = CDKSmilesGeneratorFactory.shared.newSmilesGenerator(
            flavor: [.useAromaticSymbols, .isomeric, .strict]
        )

        let smiles = generator.create(molecule)
        XCTAssertTrue(smiles.contains("."))

        let reparsed = try parser.parseSmiles(smiles)
        XCTAssertEqual(reparsed.atomCount, 2)
        XCTAssertEqual(reparsed.bondCount, 0)
    }

    func testBranchesRemainAttachedToMultivalentAtomAfterMainPath() throws {
        let molecule = try parser.parseSmiles("O=C([O-])[O-]")
        let generator = CDKSmilesGeneratorFactory.shared.newSmilesGenerator(
            flavor: [.useAromaticSymbols, .isomeric, .strict]
        )

        let smiles = generator.create(molecule)
        XCTAssertFalse(smiles.contains("C=O([O-])([O-])"))

        let reparsed = try parser.parseSmiles(smiles)
        let carbon = try XCTUnwrap(reparsed.atoms.first { $0.element.uppercased() == "C" })
        let carbonNeighbors = reparsed.neighbors(of: carbon.id)

        XCTAssertEqual(carbonNeighbors.count, 3)
        XCTAssertTrue(carbonNeighbors.allSatisfy { neighborID in
            reparsed.atom(id: neighborID)?.element.uppercased() == "O"
        })
        XCTAssertTrue(carbonNeighbors.allSatisfy { neighborID in
            reparsed.neighbors(of: neighborID).count == 1
        })
        XCTAssertEqual(reparsed.bonds(forAtom: carbon.id).filter { $0.order == .double }.count, 1)
        XCTAssertEqual(reparsed.bonds(forAtom: carbon.id).filter { $0.order == .single }.count, 2)
    }

    func testGeneratesLongLinearMoleculeWithoutRecursiveStackGrowth() throws {
        let atomCount = 5_000
        let molecule = makeLinearCarbonChain(atomCount: atomCount)
        let generator = CDKSmilesGeneratorFactory.shared.newSmilesGenerator(flavor: [.useAromaticSymbols, .strict])

        let smiles = generator.create(molecule)

        XCTAssertEqual(smiles.count, atomCount)
        XCTAssertTrue(smiles.allSatisfy { $0 == "C" })
    }

    func testGeneratesLongLinearMoleculesConcurrently() async throws {
        let molecule = makeLinearCarbonChain(atomCount: 2_000)

        try await withThrowingTaskGroup(of: String.self) { group in
            for _ in 0..<16 {
                group.addTask {
                    let generator = CDKSmilesGeneratorFactory.shared.newSmilesGenerator(
                        flavor: [.useAromaticSymbols, .strict]
                    )
                    return generator.create(molecule)
                }
            }

            for try await smiles in group {
                XCTAssertEqual(smiles.count, molecule.atomCount)
                XCTAssertTrue(smiles.allSatisfy { $0 == "C" })
            }
        }
    }

    func testNonIsomericSmilesOmitsChiralityMarker() throws {
        let molecule = try parser.parseSmiles("C[C@H](N)O")
        let generator = CDKSmilesGeneratorFactory.shared.newSmilesGenerator(flavor: [.useAromaticSymbols, .strict])

        let smiles = generator.create(molecule)
        XCTAssertFalse(smiles.contains("@"))
    }

    func testGeneratesCxAtomLabels() throws {
        let molecule = try parser.parseSmiles("CC* |$;;R1$|")
        let generator = CDKSmilesGeneratorFactory.shared.newSmilesGenerator(
            flavor: [.useAromaticSymbols, .strict, .cxAtomLabel]
        )

        XCTAssertEqual(generator.create(molecule), "CC* |$;;R1$|")
    }

    func testGeneratesCxAtomValues() throws {
        let molecule = try parser.parseSmiles("N1CN=CC1 |$_AV:HydDonor;;HydAcceptor$|")
        let generator = CDKSmilesGeneratorFactory.shared.newSmilesGenerator(
            flavor: [.useAromaticSymbols, .strict, .cxAtomValue]
        )

        let smiles = generator.create(molecule)
        XCTAssertEqual(smiles, "N1CN=CC1 |$_AV:HydDonor;;HydAcceptor$|")
        let reparsed = try parser.parseSmiles(smiles)
        XCTAssertEqual(reparsed.atoms[0].atomValue, "HydDonor")
        XCTAssertEqual(reparsed.atoms[2].atomValue, "HydAcceptor")
    }

    func testGeneratesCxCoordinates() throws {
        let molecule = try parser.parseSmiles("CCO |(,,;1,1,;2,2,)|")
        let generator = CDKSmilesGeneratorFactory.shared.newSmilesGenerator(
            flavor: [.useAromaticSymbols, .strict, .cxCoordinates]
        )

        XCTAssertEqual(generator.create(molecule), "CCO |(,,;1,1,;2,2,)|")
    }

    func testGeneratesCxMulticenter() throws {
        let molecule = try parser.parseSmiles("c1ccccc1.*Cl |m:6:0.1.2.3.4.5|")
        let generator = CDKSmilesGeneratorFactory.shared.newSmilesGenerator(
            flavor: [.useAromaticSymbols, .strict, .cxMulticenter]
        )

        XCTAssertEqual(generator.create(molecule), "c1ccccc1.*Cl |m:6:0.1.2.3.4.5|")
    }

    func testGeneratesCxPolymerSgroups() throws {
        let molecule = try parser.parseSmiles("CCCOCCO |Sg:n:1,2,3::ht|")
        let generator = CDKSmilesGeneratorFactory.shared.newSmilesGenerator(
            flavor: [.useAromaticSymbols, .strict, .cxPolymer]
        )

        XCTAssertEqual(generator.create(molecule), "CCCOCCO |Sg:n:1,2,3:n:ht|")
    }

    func testGeneratesCxRadicals() throws {
        let molecule = try parser.parseSmiles("[C]1C[CH][CH]OC1 |^1:2,3,^2:0|")
        let generator = CDKSmilesGeneratorFactory.shared.newSmilesGenerator(
            flavor: [.useAromaticSymbols, .strict, .cxRadical]
        )

        XCTAssertEqual(generator.create(molecule), "[C]1C[CH][CH]OC1 |^1:2,3,^2:0|")
    }

    func testGeneratesCxLigandOrdering() throws {
        let molecule = try parser.parseSmiles("Cl[*](Br)I |$;_R1;;$,LO:1:0.2.3|")
        let generator = CDKSmilesGeneratorFactory.shared.newSmilesGenerator(
            flavor: [.useAromaticSymbols, .strict, .cxsmiles]
        )

        let smiles = generator.create(molecule)
        XCTAssertTrue(smiles.contains("|$;R1$"))
        XCTAssertTrue(smiles.contains("LO:1:0.2.3"))

        let reparsed = try parser.parseSmiles(smiles)
        XCTAssertEqual(reparsed.atoms[1].ligandOrderingAtomIDs, [reparsed.atoms[0].id, reparsed.atoms[2].id, reparsed.atoms[3].id])
    }

    func testOmitsCoordinatesWhenFlavorDoesNotRequestCx() throws {
        let molecule = try parser.parseSmiles("CCO |(,,;1,1,;2,2,)|")
        let generator = CDKSmilesGeneratorFactory.shared.newSmilesGenerator(flavor: [.useAromaticSymbols, .strict])

        XCTAssertEqual(generator.create(molecule), "CCO")
    }

    func testOmitsCoordinatesWhenNoneArePresent() throws {
        let molecule = try parser.parseSmiles("CCO")
        let generator = CDKSmilesGeneratorFactory.shared.newSmilesGenerator(flavor: [.useAromaticSymbols, .strict, .cxCoordinates])

        XCTAssertEqual(generator.create(molecule), "CCO")
    }

    func testGeneratesReactionFragmentGroupsBeforeAtomLabels() throws {
        let reaction = try parser.parseReactionSmiles("C.*.C>> |$;R1;$,f:0.2|")
        let generator = CDKSmilesGeneratorFactory.shared.newSmilesGenerator(
            flavor: [.useAromaticSymbols, .isomeric, .strict, .cxAtomLabel, .cxFragmentGroup]
        )

        let smiles = generator.create(reaction)
        XCTAssertTrue(smiles.hasSuffix(" |f:0.2,$;R1$|"))

        let reparsed = try parser.parseReactionSmiles(smiles)
        XCTAssertEqual(reparsed.reactantCount, 2)
        XCTAssertEqual(reparsed.reactants.map(\.atomCount).sorted(), [1, 2])
    }

    func testGeneratesReactionDataSgroupsWithoutPolymerFlag() throws {
        let reaction = try parser.parseReactionSmiles("C1=CC=CC=C1C(Br)C(=O)Cl>O=C(O)c1ccccc1O>C=1C=CC=C(C1)C(C(=O)OC=2C(=CC=CC2)C(O)=O)Br |SgD::solvent:DIPEA|")
        let generator = CDKSmilesGeneratorFactory.shared.newSmilesGenerator(flavor: [.strict, .cxDataSgroups])

        let smiles = generator.create(reaction)
        XCTAssertTrue(smiles.contains("|SgD::solvent:DIPEA::|"))
    }

    func testGeneratesTersePolymerSgroup() throws {
        let molecule = try parser.parseSmiles("C1NCNC1 |Sg:n:2|")
        let generator = CDKSmilesGeneratorFactory.shared.newSmilesGenerator(flavor: .cdkDefault)

        XCTAssertEqual(generator.create(molecule), "C1NCNC1 |Sg:n:2:n:|")
    }

    func testGeneratesSimpleLinkNodeAsRepeatUnitSgroup() throws {
        let molecule = try parser.parseSmiles("C1NCNC1 |LN:2:1.3|")
        let generator = CDKSmilesGeneratorFactory.shared.newSmilesGenerator(flavor: .cdkDefault)

        XCTAssertEqual(generator.create(molecule), "C1NCNC1 |Sg:n:2:1-3:ht|")
    }

    func testGeneratesMultipleLinkNodesAsRepeatUnitSgroups() throws {
        let molecule = try parser.parseSmiles("C1NCNC1 |LN:2:1.3,3:1.4|")
        let generator = CDKSmilesGeneratorFactory.shared.newSmilesGenerator(flavor: .cdkDefault)

        XCTAssertEqual(generator.create(molecule), "C1NCNC1 |Sg:n:2:1-3:ht,Sg:n:3:1-4:ht|")
    }

    func testGeneratesLinkNodeWithAtomLabelsInCDKLayerOrder() throws {
        let molecule = try parser.parseSmiles("C1NCNC1 |LN:2:1.3,$R1$|")
        let generator = CDKSmilesGeneratorFactory.shared.newSmilesGenerator(flavor: .cdkDefault)

        XCTAssertEqual(generator.create(molecule), "*1NCNC1 |$R1$,Sg:n:2:1-3:ht|")
    }

    func testGeneratesLinkNodeWithAttachedAtomsAsExpandedRepeatUnit() throws {
        let generator = CDKSmilesGeneratorFactory.shared.newSmilesGenerator(flavor: .cdkDefault)

        let inferred = try parser.parseSmiles("C1NC(O)NC1 |LN:2:1.3|")
        let inferredSmiles = generator.create(inferred)
        XCTAssertTrue(inferredSmiles.hasSuffix(" |Sg:n:2,3:1-3:ht|"))

        let explicit = try parser.parseSmiles("C1NC(O)NC1 |LN:2:1.3.1.3|")
        let explicitSmiles = generator.create(explicit)
        XCTAssertTrue(explicitSmiles.hasSuffix(" |Sg:n:2,3:1-3:ht|"))
    }

    func testGeneratesSgroupHierarchy() throws {
        let molecule = try parser.parseSmiles("CN1CCCCC1.CO.O |Sg:c:0,1,2,3,4,5,6::,Sg:c:7,8::,Sg:c:9::,Sg:mix:0,1,2,3,4,5,6,7,8,9::,Sg:mix:7,8,9::,SgH:3:4.0,4:2.1|")
        let generator = CDKSmilesGeneratorFactory.shared.newSmilesGenerator(flavor: .cdkDefault)

        XCTAssertEqual(generator.create(molecule), "CN1CCCCC1.CO.O |Sg:c:0,1,2,3,4,5,6:c:,Sg:c:7,8:c:,Sg:c:9:c:,Sg:mix:0,1,2,3,4,5,6,7,8,9:mix:,Sg:mix:7,8,9:mix:,SgH:3:0.4,4:1.2|")
    }

    func testCollapsesSingleRacemicStereoGroupToRFlag() throws {
        let molecule = try parser.parseSmiles("C[C@H](O)[C@H](O)C1CCCCC1 |&1:1,3|")
        let generator = CDKSmilesGeneratorFactory.shared.newSmilesGenerator(
            flavor: [.useAromaticSymbols, .isomeric, .strict, .cxEnhancedStereo]
        )

        XCTAssertTrue(generator.create(molecule).hasSuffix(" |r|"))
    }

    func testCollapsesSingleRacemicStereoGroupToRFlagRegardlessOfOriginalGroupNumber() throws {
        let molecule = try parser.parseSmiles("C[C@H](O)[C@H](O)C1CCCCC1 |&3:1,3|")
        let generator = CDKSmilesGeneratorFactory.shared.newSmilesGenerator(
            flavor: [.useAromaticSymbols, .isomeric, .strict, .cxEnhancedStereo]
        )

        XCTAssertTrue(generator.create(molecule).hasSuffix(" |r|"))
    }

    func testOmitsEnhancedStereoWhenAllAbsolute() throws {
        let molecule = try parser.parseSmiles("C[C@H](O)[C@H](O)C1CCCCC1 |a:1,3|")
        let generator = CDKSmilesGeneratorFactory.shared.newSmilesGenerator(
            flavor: [.useAromaticSymbols, .isomeric, .strict, .cxEnhancedStereo]
        )

        XCTAssertFalse(generator.create(molecule).contains(" |"))
    }

    func testPreservesMixedEnhancedStereoGroups() throws {
        let molecule = try parser.parseSmiles("C[C@H](O)[C@H](O)C1CCCCC1 |o1:1,&1:3|")
        let generator = CDKSmilesGeneratorFactory.shared.newSmilesGenerator(
            flavor: [.useAromaticSymbols, .isomeric, .strict, .cxEnhancedStereo]
        )

        XCTAssertTrue(generator.create(molecule).hasSuffix(" |o1:1,&1:3|"))
    }

    func testPreservesRelativeAndAbsoluteStereoGroups() throws {
        let molecule = try parser.parseSmiles("C[C@H](O)[C@H](O)C1CCCCC1 |o1:1,a:3|")
        let generator = CDKSmilesGeneratorFactory.shared.newSmilesGenerator(
            flavor: [.useAromaticSymbols, .isomeric, .strict, .cxEnhancedStereo]
        )

        XCTAssertTrue(generator.create(molecule).hasSuffix(" |o1:1,a:3|"))
    }

    func testPromotesUngroupedStereoCentersToAbsoluteWhenMixed() throws {
        let molecule = try parser.parseSmiles("C[C@H](O)[C@H](O)C1CCCCC1 |&5:3|")
        let generator = CDKSmilesGeneratorFactory.shared.newSmilesGenerator(
            flavor: [.useAromaticSymbols, .isomeric, .strict, .cxEnhancedStereo]
        )

        XCTAssertTrue(generator.create(molecule).hasSuffix(" |a:1,&1:3|"))
    }

    private func makeLinearCarbonChain(atomCount: Int) -> Molecule {
        var molecule = Molecule(name: "Linear C\(atomCount)")
        molecule.atoms = (0..<atomCount).map { index in
            Atom(
                id: index + 1,
                element: "C",
                position: CGPoint(x: Double(index), y: 0))
        }
        molecule.bonds = (1..<atomCount).map { index in
            Bond(
                id: index,
                a1: index,
                a2: index + 1,
                order: .single)
        }
        return molecule
    }
}
