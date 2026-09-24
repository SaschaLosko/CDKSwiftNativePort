import XCTest
@testable import CDKSwiftNativePort

final class ExternalChemDrawExtensionTests: XCTestCase {
    func testChemDrawFormatsAreNotAdvertisedOrEncoded() throws {
        let molecule = try CDKSmilesParser().parseSmiles("CCO")
        for format in [CDKFileExportFormat.cdx, .cdxml] {
            XCTAssertFalse(CDKFileExporter.formats.contains { $0.format == format })
            XCTAssertNil(CDKFileExporter.format(forFileExtension: format.rawValue))
            XCTAssertThrowsError(try CDKFileExporter.writeData(molecule: molecule, as: format))
        }
        XCTAssertThrowsError(try CDKChemDrawWriter.cdx(molecules: [molecule]))
        XCTAssertThrowsError(try CDKChemDrawWriter.cdxml(molecules: [molecule]))
    }
}
