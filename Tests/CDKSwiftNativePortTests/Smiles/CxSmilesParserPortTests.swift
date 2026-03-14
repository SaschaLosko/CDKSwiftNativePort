import XCTest
@testable import CDKSwiftNativePort

final class CxSmilesParserPortTests: XCTestCase {

    func testParsesCxAtomLabelsLayer() throws {
        let parser = CDKSmilesParserFactory.shared.newSmilesParser(flavor: .cdkDefault)
        let molecule = try parser.parseSmiles("C* |$;R1$|")

        XCTAssertEqual(molecule.atomCount, 2)
        XCTAssertEqual(molecule.atoms[0].element, "C")
        XCTAssertEqual(molecule.atoms[1].element, "R1")
    }

    func testParsesCxRacemicLayerWithoutFailure() throws {
        let parser = CDKSmilesParserFactory.shared.newSmilesParser(flavor: .cdkDefault)
        let molecule = try parser.parseSmiles("CC[C@H](O)[C@H](O)CCCCCC |r|")

        XCTAssertEqual(molecule.atomCount, 12)
        XCTAssertTrue(molecule.atoms.contains(where: { $0.chirality != .none }))
    }

    func testParsesCxRacemicFragmentIndicesLayer() throws {
        let split = try CDKCxSmilesParser.split("CC[C@H](O)[C@H](O)CCCCCC |r:1|", enabled: true)
        XCTAssertEqual(split.state.racemicFragments, [1])
    }

    func testParsesCxFragmentLayerWithoutFailure() throws {
        let parser = CDKSmilesParserFactory.shared.newSmilesParser(flavor: .cdkDefault)
        let molecule = try parser.parseSmiles("C.*.C |f:0.2|")
        XCTAssertEqual(molecule.atomCount, 3)
    }

    func testParsesCxWithTitle() throws {
        let parser = CDKSmilesParserFactory.shared.newSmilesParser(flavor: .cdkDefault)
        let molecule = try parser.parseSmiles("CC |$C;C$| ethane")
        XCTAssertEqual(molecule.name, "ethane")
    }

    func testParsesCxWithoutTrailingTitleAsExplicitEmptyName() throws {
        let parser = CDKSmilesParserFactory.shared.newSmilesParser(flavor: .cdkDefault)
        let molecule = try parser.parseSmiles("CC |$C;C$|")
        XCTAssertEqual(molecule.name, "")
    }

    func testTreatsMalformedCxLayerAsTitleLikeCDK() throws {
        let parser = CDKSmilesParserFactory.shared.newSmilesParser(flavor: .cdkDefault)
        let molecule = try parser.parseSmiles("CC |$;R1$")
        XCTAssertEqual(molecule.name, "|$;R1$")
    }

    func testTreatsMalformedCxLayerBodyAsTitleLikeCDK() throws {
        let split = try CDKCxSmilesParser.split("CC |r:a|", enabled: true)
        XCTAssertEqual(split.coreSmiles, "CC")
        XCTAssertEqual(split.title, "|r:a|")
        XCTAssertTrue(split.state.atomLabels.isEmpty)
    }

    func testParsesMarkushRGroupsFromUserExample() throws {
        let parser = CDKSmilesParserFactory.shared.newSmilesParser(flavor: .cdkDefault)
        let molecule = try parser.parseSmiles("C1CNCC(*)C1 |$;;;;;R1$,LN:0:1.3,RG:_R1={OC},{Cl},{C#N}|")

        XCTAssertEqual(connectedComponents(in: molecule).count, 4)
        XCTAssertEqual(molecule.atoms.filter { $0.rGroupMembership == nil && $0.symbolToDraw == "R1" }.count, 1)
        XCTAssertEqual(molecule.atoms.filter { $0.rGroupMembership == "R1" && $0.attachmentPoint == 1 }.count, 3)
        XCTAssertEqual(molecule.sgroups.count, 1)
        XCTAssertEqual(molecule.sgroups.first?.subscriptText, "1-3")
    }

    func testParsesNestedCxInRGroupDefinitions() throws {
        let parser = CDKSmilesParserFactory.shared.newSmilesParser(flavor: .cdkDefault)
        let molecule = try parser.parseSmiles("C* |$;R1$,RG:_R1={Cl* |$;_AP1$|},{Br* |$;_AP1$|}|")

        XCTAssertEqual(connectedComponents(in: molecule).count, 3)
        XCTAssertEqual(molecule.atoms.filter { $0.rGroupMembership == "R1" && $0.attachmentPoint == 1 }.count, 2)
        XCTAssertTrue(molecule.atoms.contains { $0.rGroupMembership == "R1" && $0.element == "Cl" })
        XCTAssertTrue(molecule.atoms.contains { $0.rGroupMembership == "R1" && $0.element == "Br" })
    }

    func testParsesSimpleLinkNodeIntoRepeatUnitAnnotation() throws {
        let parser = CDKSmilesParserFactory.shared.newSmilesParser(flavor: .cdkDefault)
        let molecule = try parser.parseSmiles("C1NCNC1 |LN:2:1.3|")

        XCTAssertEqual(molecule.sgroups.count, 1)
        XCTAssertEqual(molecule.sgroups[0].kind, .structureRepeatUnit)
        XCTAssertEqual(molecule.sgroups[0].subscriptText, "1-3")
        XCTAssertEqual(molecule.sgroups[0].atomIDs.count, 1)
    }

    func testParsesMultipleLinkNodes() throws {
        let parser = CDKSmilesParserFactory.shared.newSmilesParser(flavor: .cdkDefault)
        let molecule = try parser.parseSmiles("C1NCNC1 |LN:2:1.3,3:1.4|")

        XCTAssertEqual(molecule.sgroups.count, 2)
        XCTAssertEqual(Set(molecule.sgroups.map(\.subscriptText)), Set(["1-3", "1-4"]))
    }

    func testParsesLinkNodeWithAttachedAtoms() throws {
        let parser = CDKSmilesParserFactory.shared.newSmilesParser(flavor: .cdkDefault)

        let molecule = try parser.parseSmiles("C1NC(O)NC1 |LN:2:1.3|")
        XCTAssertEqual(molecule.sgroups.count, 1)
        let linkedElements = molecule.sgroups[0].atomIDs.compactMap { atomID in
            molecule.atom(id: atomID)?.element
        }
        XCTAssertEqual(Set(linkedElements), Set(["C", "O"]))

        let withExplicitCrossingBonds = try parser.parseSmiles("C1NC(O)NC1 |LN:2:1.3.1.3|")
        XCTAssertEqual(withExplicitCrossingBonds.sgroups.count, 1)
        XCTAssertEqual(withExplicitCrossingBonds.sgroups[0].atomIDs.count, 2)
        XCTAssertEqual(Set(withExplicitCrossingBonds.sgroups[0].crossingBondIDs).count, 2)
    }

    func testPreservesDisconnectedMarkushAlternativesUsingComponentGroupIDs() throws {
        let parser = CDKSmilesParserFactory.shared.newSmilesParser(flavor: .cdkDefault)
        let molecule = try parser.parseSmiles("*c1ncccc1 |$R1$,RG:_R1={Cl},{Br},{C(=O)*.Cl |$;;R2$,RG:_R2={OC},{N}|}|")

        let r1Atoms = molecule.atoms.filter { $0.rGroupMembership == "R1" }
        let grouped = Dictionary(grouping: r1Atoms, by: { $0.componentGroupID ?? -1 })

        XCTAssertGreaterThanOrEqual(grouped.count, 3)
        XCTAssertTrue(grouped.values.contains { atoms in
            atoms.contains(where: { $0.element == "Cl" }) &&
            atoms.contains(where: { $0.element == "C" || $0.element == "O" })
        })
    }

    func testTreatsNonCxLayerAsTitleLikeCDK() throws {
        let parser = CDKSmilesParserFactory.shared.newSmilesParser(flavor: .cdkDefault)
        let molecule = try parser.parseSmiles("c1ccccc1 |<benzene>|")
        XCTAssertEqual(molecule.name, "|<benzene>|")
    }

    func testTreatsTruncatedCxPipeAsTitleLikeCDK() throws {
        let parser = CDKSmilesParserFactory.shared.newSmilesParser(flavor: .cdkDefault)
        let molecule = try parser.parseSmiles("c1ccccc1 |")
        XCTAssertEqual(molecule.name, "|")
    }

    func testParsesEscapedAtomLabelEntity() throws {
        let parser = CDKSmilesParserFactory.shared.newSmilesParser(flavor: .cdkDefault)
        let molecule = try parser.parseSmiles("**.c1ccccc1CC |$R&#39;$|")

        XCTAssertEqual(molecule.atoms[0].element, "R'")
    }

    func testParsesEscapedComplexAtomLabelEntity() throws {
        let split = try CDKCxSmilesParser.split("CCCC |$;;;&#40;C&#40;R41&#41;&#40;R41&#41;&#41;n$|", enabled: true)
        XCTAssertEqual(split.state.atomLabels[3], "(C(R41)(R41))n")
    }

    func testParsesEmptyAtomLabelLayerWithoutMutatingAtoms() throws {
        let parser = CDKSmilesParserFactory.shared.newSmilesParser(flavor: .cdkDefault)
        let molecule = try parser.parseSmiles("CC |$$|")

        XCTAssertEqual(molecule.atoms.map(\.element), ["C", "C"])
        XCTAssertEqual(molecule.name, "")
    }

    func testParsesAttachmentPointPseudoLabel() throws {
        let parser = CDKSmilesParserFactory.shared.newSmilesParser(flavor: .cdkDefault)
        let molecule = try parser.parseSmiles("**.c1ccccc1CC |$;;;;;;;;;_AP1$|")

        XCTAssertEqual(molecule.atoms[9].attachmentPoint, 1)
        XCTAssertEqual(molecule.atoms[9].explicitHydrogenCount, 0)
    }

    func testParsesCxAtomValuesLayer() throws {
        let parser = CDKSmilesParserFactory.shared.newSmilesParser(flavor: .cdkDefault)
        let molecule = try parser.parseSmiles("N1CN=CC1 |$_AV:HydDonor;;HydAcceptor$|")

        XCTAssertEqual(molecule.atoms[0].atomValue, "HydDonor")
        XCTAssertEqual(molecule.atoms[2].atomValue, "HydAcceptor")
    }

    func testParses2DCoordinatesIntoAtoms() throws {
        let parser = CDKSmilesParserFactory.shared.newSmilesParser(flavor: .cdkDefault)
        let molecule = try parser.parseSmiles("CCC |(0,1,;0,2,;0,3,)|")

        XCTAssertEqual(molecule.atoms[0].position.x, 0, accuracy: 0.0001)
        XCTAssertEqual(molecule.atoms[0].position.y, 1, accuracy: 0.0001)
        XCTAssertNil(molecule.atoms[0].zPosition)
        XCTAssertEqual(molecule.atoms[2].position.y, 3, accuracy: 0.0001)
    }

    func testParses3DCoordinatesIntoAtoms() throws {
        let parser = CDKSmilesParserFactory.shared.newSmilesParser(flavor: .cdkDefault)
        let molecule = try parser.parseSmiles("CCC |(0,1,1;0,2,1;0,3,1)|")

        XCTAssertEqual(molecule.atoms[0].position.y, 1, accuracy: 0.0001)
        XCTAssertEqual(try XCTUnwrap(molecule.atoms[0].zPosition), 1, accuracy: 0.0001)
        XCTAssertEqual(try XCTUnwrap(molecule.atoms[2].zPosition), 1, accuracy: 0.0001)
    }

    func testParsesPositionalVariationIntoExtMulticenterSgroup() throws {
        let parser = CDKSmilesParserFactory.shared.newSmilesParser(flavor: .cdkDefault)
        let molecule = try parser.parseSmiles("**.c1ccccc1CC |m:1:2.3.4.5.6.7|")

        let sgroup = try XCTUnwrap(molecule.sgroups.first(where: { $0.kind == .extMulticenter }))
        XCTAssertEqual(sgroup.keyword, "m")
        XCTAssertEqual(sgroup.atomIDs.count, 7)
        XCTAssertEqual(sgroup.crossingBondIDs.count, 1)
    }

    func testParsesPolymerSgroup() throws {
        let parser = CDKSmilesParserFactory.shared.newSmilesParser(flavor: .cdkDefault)
        let molecule = try parser.parseSmiles("**.c1ccccc1CC |Sg:n:8:m:ht|")

        let sgroup = try XCTUnwrap(molecule.sgroups.first)
        XCTAssertEqual(sgroup.kind, .structureRepeatUnit)
        XCTAssertEqual(sgroup.keyword, "n")
        XCTAssertEqual(sgroup.subscriptText, "m")
        XCTAssertEqual(sgroup.connectivity, "ht")
        XCTAssertEqual(sgroup.atomIDs.count, 1)
        XCTAssertEqual(sgroup.crossingBondIDs.count, 2)
    }

    func testParsesTersePolymerSgroup() throws {
        let parser = CDKSmilesParserFactory.shared.newSmilesParser(flavor: .cdkDefault)
        let molecule = try parser.parseSmiles("C1NCNC1 |Sg:n:2|")

        let sgroup = try XCTUnwrap(molecule.sgroups.first)
        XCTAssertEqual(sgroup.keyword, "n")
        XCTAssertEqual(sgroup.subscriptText, "n")
        XCTAssertEqual(sgroup.connectivity, "")
        XCTAssertEqual(sgroup.atomIDs.count, 1)
    }

    func testParsesDataSgroup() throws {
        let parser = CDKSmilesParserFactory.shared.newSmilesParser(flavor: .cdkDefault)
        let molecule = try parser.parseSmiles("C1=CC=CC=C1C(Br)C(=O)Cl |SgD::solvent:DIPEA::|")

        let sgroup = try XCTUnwrap(molecule.sgroups.first(where: { $0.kind == .data }))
        XCTAssertEqual(sgroup.dataFieldName, "solvent")
        XCTAssertEqual(sgroup.dataValue, "DIPEA")
        XCTAssertEqual(sgroup.dataOperator, "")
        XCTAssertEqual(sgroup.dataUnit, "")
    }

    func testParsesDataSgroupWithEscapedCDKFieldAlias() throws {
        let parser = CDKSmilesParserFactory.shared.newSmilesParser(flavor: .cdkDefault)
        let molecule = try parser.parseSmiles("C |SgD::cdk&#58;Arrow:RES|")

        let sgroup = try XCTUnwrap(molecule.sgroups.first(where: { $0.kind == .data }))
        XCTAssertEqual(sgroup.dataFieldName, "cdk:Arrow")
        XCTAssertEqual(sgroup.dataValue, "RES")
    }

    func testParsesDataSgroupWithDotSeparatedCDKFieldAlias() throws {
        let parser = CDKSmilesParserFactory.shared.newSmilesParser(flavor: .cdkDefault)
        let molecule = try parser.parseSmiles("C |SgD::cdk.Arrow:RES|")

        let sgroup = try XCTUnwrap(molecule.sgroups.first(where: { $0.kind == .data }))
        XCTAssertEqual(sgroup.dataFieldName, "cdk:Arrow")
        XCTAssertEqual(sgroup.dataValue, "RES")
    }

    func testParsesSgroupHierarchy() throws {
        let parser = CDKSmilesParserFactory.shared.newSmilesParser(flavor: .cdkDefault)
        let molecule = try parser.parseSmiles("CN1CCCCC1.CO.O |Sg:c:0,1,2,3,4,5,6::,Sg:c:7,8::,Sg:c:9::,Sg:mix:0,1,2,3,4,5,6,7,8,9::,Sg:mix:7,8,9::,SgH:3:4.0,4:2.1|")

        XCTAssertEqual(molecule.sgroups.count, 5)
        XCTAssertEqual(Set(molecule.sgroups[3].childGroupIndices), Set([0, 4]))
        XCTAssertEqual(Set(molecule.sgroups[4].childGroupIndices), Set([1, 2]))
    }

    func testParsesPositionalVariationImpliedLayerEntries() throws {
        let split = try CDKCxSmilesParser.split("CC |m:2:5.6.7.8.9.10,4:5.6.7.8.9|", enabled: true)

        XCTAssertEqual(split.state.positionalVariations[2], [5, 6, 7, 8, 9, 10])
        XCTAssertEqual(split.state.positionalVariations[4], [5, 6, 7, 8, 9])
    }

    func testParsesStereoGroupState() throws {
        let split = try CDKCxSmilesParser.split("CC[C@H](O)[C@H](O)CCCCCC |o1:2,&5:4,a:1|", enabled: true)

        XCTAssertEqual(split.state.stereoGroups[2], CDKCxSmilesParser.encodeStereoGroup(kind: "or", number: 1))
        XCTAssertEqual(split.state.stereoGroups[4], CDKCxSmilesParser.encodeStereoGroup(kind: "and", number: 5))
        XCTAssertEqual(split.state.stereoGroups[1], CDKCxSmilesParser.encodeStereoGroup(kind: "abs", number: 0))
    }

    func testParsesSingleAndStereoGroupState() throws {
        let split = try CDKCxSmilesParser.split("CC |&1:0,1|", enabled: true)
        XCTAssertEqual(split.state.stereoGroups[0], CDKCxSmilesParser.encodeStereoGroup(kind: "and", number: 1))
        XCTAssertEqual(split.state.stereoGroups[1], CDKCxSmilesParser.encodeStereoGroup(kind: "and", number: 1))
    }

    func testParsesSingleOrStereoGroupState() throws {
        let split = try CDKCxSmilesParser.split("CC |o1:0,1|", enabled: true)
        XCTAssertEqual(split.state.stereoGroups[0], CDKCxSmilesParser.encodeStereoGroup(kind: "or", number: 1))
        XCTAssertEqual(split.state.stereoGroups[1], CDKCxSmilesParser.encodeStereoGroup(kind: "or", number: 1))
    }

    func testParsesRadicalLayerIntoAtoms() throws {
        let parser = CDKSmilesParserFactory.shared.newSmilesParser(flavor: .cdkDefault)
        let molecule = try parser.parseSmiles("[N]1C=CC=C1 |c:1,3,^1:0|")

        XCTAssertEqual(molecule.atoms[0].radicalType, .monovalent)
        XCTAssertEqual(molecule.atoms[0].radical, 1)
    }

    func testParsesDivalentRadicalLayerIntoAtoms() throws {
        let parser = CDKSmilesParserFactory.shared.newSmilesParser(flavor: .cdkDefault)
        let molecule = try parser.parseSmiles("[C]1C2=CC=CC=C2C2=CC=CC=C12 |c:3,5,10,t:1,8,12,^3:0|")

        XCTAssertEqual(molecule.atoms[0].radicalType, .divalentSinglet)
        XCTAssertEqual(molecule.atoms[0].radical, 2)
    }

    func testParsesLigandOrderingIntoAnchorAtom() throws {
        let parser = CDKSmilesParserFactory.shared.newSmilesParser(flavor: .cdkDefault)
        let molecule = try parser.parseSmiles("Cl[*](Br)I |$;_R1;;$,LO:1:0.2.3|")

        XCTAssertEqual(molecule.atoms[1].element, "R1")
        XCTAssertEqual(molecule.atoms[1].ligandOrderingAtomIDs, [molecule.atoms[0].id, molecule.atoms[2].id, molecule.atoms[3].id])
    }

    func testParsesHighlightAndWedgeLayersWithoutFailure() throws {
        let split = try CDKCxSmilesParser.split("C1NCNC1 |ha:0,1,3,hb:2,4,wD:1.0,wU:2.1|", enabled: true)

        XCTAssertEqual(split.state.atomHighlights, [0, 1, 3])
        XCTAssertEqual(split.state.bondHighlights, [2, 4])
        XCTAssertEqual(split.state.bondDisplays.count, 2)
    }

    func testParsesHighlightedAtomsAndBondsIntoMoleculeSelectionState() throws {
        let parser = CDKSmilesParserFactory.shared.newSmilesParser(flavor: .cdkDefault)
        let molecule = try parser.parseSmiles("C1NCNC1 |ha:0,1,3,hb:2,4|")

        XCTAssertEqual(molecule.highlightedAtomIDs,
                       [molecule.atoms[0].id, molecule.atoms[1].id, molecule.atoms[3].id])
        XCTAssertEqual(molecule.highlightedBondIDs,
                       [molecule.bonds[2].id, molecule.bonds[4].id])
    }

    func testParsesAtomLabelsAndAtomValuesTogether() throws {
        let parser = CDKSmilesParserFactory.shared.newSmilesParser(flavor: .cdkDefault)
        let molecule = try parser.parseSmiles("[H]C1=C([H])N2C(=O)C(=C([O-])[N+](CC3=CN=C(Cl)S3)=C2C(C)=C1[H])C1=CC(*)=CC=C1.** |$;;;;;;;;;;;;;;;;;;;;;;;;;;R;;;;RA;$,$_AV:;;;;;;;;;;;;;;;;;;;;;;;;2;;;4;5;6;;$,m:31:29.28.27.25.24.23|")

        XCTAssertEqual(molecule.atoms[26].element, "R")
        XCTAssertEqual(molecule.atoms[30].element, "RA")
        XCTAssertEqual(molecule.atoms[24].atomValue, "2")
        XCTAssertEqual(molecule.atoms[27].atomValue, "4")
        XCTAssertEqual(molecule.atoms[28].atomValue, "5")
        XCTAssertEqual(molecule.atoms[29].atomValue, "6")
        XCTAssertEqual(molecule.sgroups.filter { $0.kind == .extMulticenter }.count, 1)
    }

    func testPreservesVariableAttachCrossingBondsBugfix() throws {
        let parser = CDKSmilesParserFactory.shared.newSmilesParser(flavor: .cdkDefault)
        let molecule = try parser.parseSmiles("[H]OCCO.C* |lp:1:2,4:2,m:6:3.2,Sg:n:1,2,3,5::ht|")

        let sgroup = try XCTUnwrap(molecule.sgroups.first(where: { $0.kind == .structureRepeatUnit }))
        XCTAssertEqual(sgroup.atomIDs.count, 4)
        XCTAssertEqual(sgroup.crossingBondIDs.count, 2)
    }

    private func connectedComponents(in molecule: Molecule) -> [Set<Int>] {
        var adjacency: [Int: Set<Int>] = [:]
        for atom in molecule.atoms {
            _ = adjacency[atom.id, default: []]
        }
        for bond in molecule.bonds {
            adjacency[bond.a1, default: []].insert(bond.a2)
            adjacency[bond.a2, default: []].insert(bond.a1)
        }

        var seen = Set<Int>()
        var components: [Set<Int>] = []

        for atom in molecule.atoms where !seen.contains(atom.id) {
            var queue = [atom.id]
            var component = Set<Int>()
            seen.insert(atom.id)

            while let current = queue.popLast() {
                component.insert(current)
                for neighbor in adjacency[current, default: []] where !seen.contains(neighbor) {
                    seen.insert(neighbor)
                    queue.append(neighbor)
                }
            }

            components.append(component)
        }

        return components
    }
}
