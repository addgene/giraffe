load "features/urls.rb"

# Commonly used webrat steps
# http://github.com/brynary/webrat

Given /^that I am on (.+)$/ do |page_name|
  visit $urls[page_name]
end

When /^I goto (.+)$/ do |page_name|
  visit $urls[page_name]
end

When /^I visit "([^\"]*)"$/ do |url|
  visit url
end

When /^I click on "([^\"]*)"$/ do |link|
  click_link(link)
end

When /^I click on link "([^\"]*)"$/ do |link|
  click_link(link)
end

When /^I press "([^\"]*)"$/ do |button|
  click_button(button)
end

When /^I follow "([^\"]*)"$/ do |link|
  click_link(link)
end

When /^I fill in "([^\"]*)" with random number$/ do |field|
  value = (0...5).map{48.+(rand(10)).chr}.join
  fill_in(field, :with => value) 
end

When /^I fill in "([^\"]*)" with random string$/ do |field|
  value = (0...16).map{97.+(rand(26)).chr}.join
  fill_in(field, :with => value) 
end

When /^I fill in "([^\"]*)" with "([^\"]*)"$/ do |field, value|
  fill_in(field, :with => value) 
end

When /^I select "([^\"]*)" from "([^\"]*)"$/ do |value, field|
  select(value, :from => field) 
end

# Use this step in conjunction with Rail's datetime_select helper. For example:
# When I select "December 25, 2008 10:00" as the date and time 
When /^I select "([^\"]*)" as the date and time$/ do |time|
  select_datetime(time)
end

# Use this step when using multiple datetime_select helpers on a page or 
# you want to specify which datetime to select. Given the following view:
#   <%%= f.label :preferred %><br />
#   <%%= f.datetime_select :preferred %>
#   <%%= f.label :alternative %><br />
#   <%%= f.datetime_select :alternative %>
# The following steps would fill out the form:
# When I select "November 23, 2004 11:20" as the "Preferred" data and time
# And I select "November 25, 2004 10:30" as the "Alternative" data and time
When /^I select "([^\"]*)" as the "([^\"]*)" date and time$/ do |datetime, datetime_label|
  select_datetime(datetime, :from => datetime_label)
end

# Use this step in conjunction with Rail's time_select helper. For example:
# When I select "2:20PM" as the time
# Note: Rail's default time helper provides 24-hour time-- not 12 hour time. Webrat
# will convert the 2:20PM to 14:20 and then select it. 
When /^I select "([^\"]*)" as the time$/ do |time|
  select_time(time)
end

# Use this step when using multiple time_select helpers on a page or you want to
# specify the name of the time on the form.  For example:
# When I select "7:30AM" as the "Gym" time
When /^I select "([^\"]*)" as the "([^\"]*)" time$/ do |time, time_label|
  select_time(time, :from => time_label)
end

# Use this step in conjunction with Rail's date_select helper.  For example:
# When I select "February 20, 1981" as the date
When /^I select "([^\"]*)" as the date$/ do |date|
  select_date(date)
end

# Use this step when using multiple date_select helpers on one page or
# you want to specify the name of the date on the form. For example:
# When I select "April 26, 1982" as the "Date of Birth" date
When /^I select "([^\"]*)" as the "([^\"]*)" date$/ do |date, date_label|
  select_date(date, :from => date_label)
end

When /^I check "([^\"]*)"$/ do |field|
  check(field) 
end

When /^I uncheck "([^\"]*)"$/ do |field|
  uncheck(field) 
end

When /^I choose "([^\"]*)"$/ do |field|
  choose(field)
end

When /^I attach the file at "([^\"]*)" to "([^\"]*)"$/ do |path, field|
  attach_file(field, path)
end

Then /^I should see "([^\"]*)"$/ do |text|
  response_body.should contain(text)
end

Then /^I should not see "([^\"]*)"$/ do |text|
  response_body.should_not contain(text)
end

Then /^the "([^\"]*)" field should contain "([^\"]*)"$/ do |field, value|
  field_labeled(field).value.should =~ /#{value}/
end

Then /^the "([^\"]*)" field should not contain "([^\"]*)"$/ do |field, value|
  field_labeled(field).value.should_not =~ /#{value}/
end
    
Then /^the "([^\"]*)" checkbox should be checked$/ do |label|
  field_labeled(label).should be_checked
end

Then /^I should be on (.+)$/ do |page_name|
  expected_path = URI.parse($urls[page_name]).select(:path, :query).compact.join('?')
  current_path = URI.parse(current_url).select(:path, :query).compact.join('?')
  current_path.should == expected_path
end

Then /^the page should have "([^\"]*)" area$/ do |area|
  response_body.should have_selector("#"+area)
end

Then /^the page should not have "([^\"]*)" area$/ do |area|
  response_body.should_not have_selector("#"+area)
end

Then /^the "([^\"]*)" area should contain "([^\"]*)"$/ do |area, value|
  response_body.should have_selector("#"+area) do |a|
    a.should contain(value)
  end
end

Then /^the "([^\"]*)" area should not contain "([^\"]*)"$/ do |area, value|
  response_body.should have_selector("#"+area) do |a|
    a.should_not contain(value)
  end
end

Then /^the "([^\"]*)" area should match "([^\"]*)"$/ do |area, value|
  response_body.should have_selector("#"+area) do |a|
    a.should == value
  end
end

