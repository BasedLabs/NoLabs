describe('Main page', () => {
  const inputAminoAcidSequenceId = '#inputSequence';
  const submitAminoAcidSequenceId = '#submitInference';
  const localisationImageClass = '.localisation-image';

  it('Opened main page', () => {
    cy.visit('http://127.0.0.1:5000');
  });

  it('Amino-acid sequence input is visible', () => {
    cy.visit('http://127.0.0.1:5000');
    cy.get(inputAminoAcidSequenceId).should('exist');
  });

  it('Inference page is opening', () => {
    cy.visit('http://127.0.0.1:5000');
    cy.get(inputAminoAcidSequenceId).click();
    cy.get(inputAminoAcidSequenceId).type(`MNGTEGPNFYVPFSNATGVVRSPFEYPQYYLAEPWQFSMLAAYMFLLIVLGFPINFLTLYVTVQHKKLRTPLNYILLNLAVADLFMVLGGFTSTLYTSLHGYFVFGPTGCNLEGFFATLGGEIALWSLVVLAIERYVVVCKPMSNFRFGENHAIMGVAFTWVMALACAAPPLAGWSRYIPEGLQCSCGIDYYTLKPEVNNESFVIYMFVVHFTIPMIIIFFCYGQLVFTVKEAAAQQQESATTQKAEKEVTRMVIIMVIAFLICWVPYASVAFYIFTHQGSNFGPIFMTIPAFFAKSAAIYNPVIYIMMNKQFRNCMLTTICCGKNPLGDDEASATVSKTETSQVAPA`);
    cy.get(submitAminoAcidSequenceId).click();
    cy.location('pathname', {timeout: 120000}).should('include', '/inference');
    cy.get(localisationImageClass).should('exist');
  });
})