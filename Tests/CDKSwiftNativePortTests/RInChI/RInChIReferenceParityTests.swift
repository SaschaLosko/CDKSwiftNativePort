import Foundation
import XCTest
@testable import CDKSwiftNativePort

final class RInChIReferenceParityTests: XCTestCase {
    private struct ReferenceCase {
        let id: String
        let inputFile: String
        let outputFile: String
        let forceEquilibrium: Bool
    }

    private struct Gap: Equatable {
        let caseID: String
        let field: String
        let expected: String
        let actual: String
    }

    private struct Report: Encodable {
        struct ReportGap: Encodable {
            let caseID: String
            let field: String
            let expected: String
            let actual: String
        }

        let suite: String
        let caseCount: Int
        let gaps: [ReportGap]
    }

    private static let cases: [ReferenceCase] = [
        .init(id: "r01", inputFile: "r01_1_struct_reactant_1_nostruct_product.rxn", outputFile: "r01_1_struct_reactant_1_nostruct_product-rinchi.txt", forceEquilibrium: false),
        .init(id: "r02", inputFile: "r02_1_struct_reactant_0_product.rxn", outputFile: "r02_1_struct_reactant_0_product-rinchi.txt", forceEquilibrium: false),
        .init(id: "r03", inputFile: "r03_1_struct_reactant_1_nostruct_product.rxn", outputFile: "r03_1_struct_reactant_1_nostruct_product-rinchi.txt", forceEquilibrium: false),
        .init(id: "r04", inputFile: "r04_1_struct_reactant_1_nostruct_product_test.rxn", outputFile: "r04_1_struct_reactant_1_nostruct_product_test-rinchi.txt", forceEquilibrium: false),
        .init(id: "r05", inputFile: "r05_0_reactant_1_struct_product.rxn", outputFile: "r05_0_reactant_1_struct_product-rinchi.txt", forceEquilibrium: false),
        .init(id: "r10", inputFile: "r10_1_struct_reactant_1_struct_product.rxn", outputFile: "r10_1_struct_reactant_1_struct_product-rinchi.txt", forceEquilibrium: false),
        .init(id: "r11", inputFile: "r11_2_struct_reactant_0_product.rxn", outputFile: "r11_2_struct_reactant_0_product-rinchi.txt", forceEquilibrium: false),
        .init(id: "r12", inputFile: "r12_1_nostruct_reactant_1_nostruct_product.rxn", outputFile: "r12_1_nostruct_reactant_1_nostruct_product-rinchi.txt", forceEquilibrium: false),
        .init(id: "r13", inputFile: "r13_1_nostruct_reactant_1_nostruct_product.rxn", outputFile: "r13_1_nostruct_reactant_1_nostruct_product-rinchi.txt", forceEquilibrium: false),
        .init(id: "r15", inputFile: "r15_4_struct_reactant_1_struct_product.rxn", outputFile: "r15_4_struct_reactant_1_struct_product-rinchi.txt", forceEquilibrium: false),
        .init(id: "r17", inputFile: "r17_rinchi_repo_5_variations_1_step_each.rxn", outputFile: "r17_rinchi_repo_5_variations_1_step_each-rinchi.txt", forceEquilibrium: false),
        .init(id: "r20", inputFile: "r20_rinchi_repo_example_03_metab_udm.rxn", outputFile: "r20_rinchi_repo_example_03_metab_udm-rinchi.txt", forceEquilibrium: false),
        .init(id: "r23", inputFile: "r23_rinchi_repo_1_nostruct_reactant_1_struct_product.rxn", outputFile: "r23_rinchi_repo_1_nostruct_reactant_1_struct_product-rinchi.txt", forceEquilibrium: false),
        .init(id: "r24", inputFile: "r24_rinchi_repo_1_struct_reactant_R_1_struct_product_X.rxn", outputFile: "r24_rinchi_repo_1_struct_reactant_R_1_struct_product_X-rinchi.txt", forceEquilibrium: false),
        .init(id: "r25", inputFile: "r25_rinchi_repo_2_reactant_asterisk_1_nostruct_product.rxn", outputFile: "r25_rinchi_repo_2_reactant_asterisk_1_nostruct_product-rinchi.txt", forceEquilibrium: false),
        .init(id: "r26", inputFile: "r26_1_nostruct_reactant_1_nostruct_product_1_nostruct_agent_test.rxn", outputFile: "r26_1_nostruct_reactant_1_nostruct_product_1_nostruct_agent_test.txt", forceEquilibrium: false),
        .init(id: "r27", inputFile: "r27_uspto_1976_US03930949_16_test.rxn", outputFile: "r27_uspto_1976_US03930949_16_test.txt", forceEquilibrium: false),
        .init(id: "r28", inputFile: "r28.rxn", outputFile: "r28.txt", forceEquilibrium: false),
        .init(id: "r29", inputFile: "r29_forceEquilibrium.rxn", outputFile: "r29_forceEquilibrium.txt", forceEquilibrium: true),
        .init(id: "r30", inputFile: "r30_malate_oxaloacetate_notChiral_forceEquilibrium.rxn", outputFile: "r30_malate_oxaloacetate_notChiral_forceEquilibrium.txt", forceEquilibrium: true),
        .init(id: "r31", inputFile: "r31_malate_oxaloacetate_chiral_forceEquilibrium.rxn", outputFile: "r31_malate_oxaloacetate_chiral_forceEquilibrium.txt", forceEquilibrium: true),
    ]

    func testReferenceCorpusMetadataIsPresent() {
        XCTAssertEqual(Self.cases.count, 21)
        XCTAssertTrue(FileManager.default.fileExists(atPath: referenceDirectory().path))
    }

    func testVendoredReferenceParity() throws {
        var gaps: [Gap] = []

        for caseEntry in Self.cases {
            let reactionText = try String(contentsOf: referenceDirectory().appendingPathComponent(caseEntry.inputFile), encoding: .utf8)
            let expected = try loadReferenceOutput(named: caseEntry.outputFile)
            let reaction = try CDKRXNReader.readReaction(text: reactionText)
            let generator = CDKRInChIGenerator(
                options: CDKRInChIOptions(forceEquilibrium: caseEntry.forceEquilibrium)
            ).generate(reaction)

            gaps.append(contentsOf: compare(caseID: caseEntry.id,
                                            field: "RInChI",
                                            expected: expected["RInChI"] ?? "",
                                            actual: generator.getRInChI() ?? ""))
            gaps.append(contentsOf: compare(caseID: caseEntry.id,
                                            field: "RAuxInfo",
                                            expected: expected["RAuxInfo"] ?? "",
                                            actual: generator.getAuxInfo() ?? ""))
            gaps.append(contentsOf: compare(caseID: caseEntry.id,
                                            field: "Long-RInChIKey",
                                            expected: expected["Long-RInChIKey"] ?? "",
                                            actual: generator.getLongRInChIKey() ?? ""))
            gaps.append(contentsOf: compare(caseID: caseEntry.id,
                                            field: "Short-RInChIKey",
                                            expected: expected["Short-RInChIKey"] ?? "",
                                            actual: generator.getShortRInChIKey() ?? ""))
            gaps.append(contentsOf: compare(caseID: caseEntry.id,
                                            field: "Web-RInChIKey",
                                            expected: expected["Web-RInChIKey"] ?? "",
                                            actual: generator.getWebRInChIKey() ?? ""))
        }

        let report = Report(
            suite: "CDK RInChI Reference Corpus",
            caseCount: Self.cases.count,
            gaps: gaps.map { .init(caseID: $0.caseID, field: $0.field, expected: $0.expected, actual: $0.actual) }
        )
        let reportPath = try writeReport(report)

        XCTAssertTrue(
            gaps.isEmpty,
            """
            CDK RInChI reference parity is not complete.
            Detailed report: \(reportPath.path)
            """
        )
    }

    private func compare(caseID: String, field: String, expected: String, actual: String) -> [Gap] {
        expected == actual ? [] : [Gap(caseID: caseID, field: field, expected: expected, actual: actual)]
    }

    private func referenceDirectory() -> URL {
        URL(fileURLWithPath: #filePath)
            .deletingLastPathComponent()
            .appendingPathComponent("UpstreamReference", isDirectory: true)
    }

    private func loadReferenceOutput(named fileName: String) throws -> [String: String] {
        let text = try String(contentsOf: referenceDirectory().appendingPathComponent(fileName), encoding: .utf8)
        var outputs: [String: String] = [:]
        for line in text
            .replacingOccurrences(of: "\r\n", with: "\n")
            .replacingOccurrences(of: "\r", with: "\n")
            .split(separator: "\n", omittingEmptySubsequences: true)
            .map(String.init) {
            if line.hasPrefix("RInChI=") {
                outputs["RInChI"] = line
            } else if line.hasPrefix("RAuxInfo=") {
                outputs["RAuxInfo"] = line
            } else if line.hasPrefix("Long-RInChIKey=") {
                outputs["Long-RInChIKey"] = line
            } else if line.hasPrefix("Short-RInChIKey=") {
                outputs["Short-RInChIKey"] = line
            } else if line.hasPrefix("Web-RInChIKey=") {
                outputs["Web-RInChIKey"] = line
            }
        }
        return outputs
    }

    private func writeReport(_ report: Report) throws -> URL {
        let encoder = JSONEncoder()
        encoder.outputFormatting = [.prettyPrinted, .sortedKeys]
        let data = try encoder.encode(report)
        let url = URL(fileURLWithPath: NSTemporaryDirectory())
            .appendingPathComponent("cdk_rinchi_reference_parity_report.json")
        try data.write(to: url, options: .atomic)
        return url
    }
}
