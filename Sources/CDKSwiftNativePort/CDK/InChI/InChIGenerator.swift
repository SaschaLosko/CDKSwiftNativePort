import Foundation

/// Swift counterpart of CDK's InChI generator facade for structure-to-InChI conversion.
/// Uses the native Swift CDK port implementation.
public final class CDKInChIGenerator {
    private let inchi: String?
    private let inchiKey: String?

    private(set) var status: CDKInChIStatus = .success
    private(set) var message: String = ""

    public init(molecule: Molecule) {
        do {
            let result = try CDKInChINativeGenerator.generate(for: molecule)
            self.inchi = result.inchi
            self.inchiKey = result.inchiKey
            self.status = result.status
            self.message = result.message
        } catch {
            self.inchi = nil
            self.inchiKey = nil
            self.status = .error
            self.message = (error as? LocalizedError)?.errorDescription ?? error.localizedDescription
        }
    }

    public func getStatus() -> CDKInChIStatus {
        status
    }

    public func getMessage() -> String {
        message
    }

    public func getInchi() throws -> String {
        guard let inchi else {
            throw ChemError.parseFailed(message.isEmpty ? "Could not generate InChI." : message)
        }
        return inchi
    }

    public func getInchiKey() throws -> String {
        guard let inchiKey else {
            throw ChemError.parseFailed(message.isEmpty ? "Could not generate InChIKey." : message)
        }
        return inchiKey
    }
}
